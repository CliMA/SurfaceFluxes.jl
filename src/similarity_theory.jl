### REVAMP 
using Statistics
import Base: -
import Statistics:norm

"""
    SimilarityScales

Holds the Monin-Obukhov solution state u★, θ★, q★
"""
struct SimilarityScales{U, T, Q}
    momentum::U # u★
    temperature::T # θ★
    water_vapor::Q # q★
end

"""
    -(a::SimilarityScales, b::SimilarityScales)

Component-wise subtraction
"""
function -(a::SimilarityScales, b::SimilarityScales)
    Δu = a.momentum - b.momentum
    Δθ = a.temperature - b.temperature
    Δq = a.water_vapor - b.water_vapor
    return SimilarityScales(Δu, Δθ, Δq)
end

Statistics.norm(a::SimilarityScales) = norm(a.momentum) + norm(a.temperature) + norm(a.water_vapor)

abstract type AbstractStateVariables end

"""
    SurfaceState

Holds the functions required to compute roughness lengths, surface conditions
Use parametric types
"""
struct SurfaceState <: AbstractStateVariables
    roughness_lengths
    u_s 
    v_s
    q_s
    θ_s
    h_s
    args
end
struct AtmosState <: AbstractStateVariables
    u_a 
    v_a
    q_a
    θ_a
    h_a
    gustiness_parameter # Uᵍᵢ
    h_boundary_layer # TODO computation
    args
end

function surface_variable(surface_state_var::Function, surface_args, similarity_scales, atmos_state, params)
    return surface_state_var(surface_args, similarity_scales, atmos_state, params)
end

function surface_variable(surface_state_var::FT, _...) where {FT <: AbstractFloat}
    return surface_state_var
end

@inline function buoyancy_scale(θ★, q★, thermo_state, thermo_params, g)
    𝒯ₐ = TD.virtual_temperature(thermo_params, thermo_state)
    qₐ = TD.vapor_specific_humidity(thermo_params, thermo_state)
    ε  = TD.Parameters.molmass_ratio(thermo_params)
    δ  = ε - 1 # typically equal to 0.608
    b★ = g / 𝒯ₐ * (θ★ * (1 + δ * qₐ) + δ * 𝒯ₐ * q★)
    return b★
end

"""
    state_differences(surface_state, atmos_state, similarity_scales, params)

returns `NamedTuple` of state differences across first interior cell
"""
function state_differences(surface_state, atmos_state, similarity_scales, params)
    surface_args = surface_state.args
    
    q_a=atmos_state.q_a
    q_s=surface_variable(surface_state.q_s, 
                         surface_args, 
                         similarity_scales, atmos_state, params)
    Δq=q_a - q_s
    
    u_a=atmos_state.u_a
    u_s=surface_variable(surface_state.u_s, 
                         surface_args, 
                         similarity_scales, atmos_state, params)
    Δu= u_a - u_s
    
    v_a=atmos_state.v_a
    v_s=surface_variable(surface_state.v_s, 
                         surface_args, 
                         similarity_scales, atmos_state, params)
    Δv=v_a - v_s
    
    θ_a=atmos_state.θ_a
    θ_s=surface_variable(surface_state.θ_s, 
                         surface_args, 
                         similarity_scales, atmos_state, params)
    Δθ=θ_a - θ_s
    
    h_a=atmos_state.h_a
    h_s=surface_state.h_s
    Δh=h_a-h_s 

    return (; Δu=Δu, Δv=Δv, Δθ=Δθ, Δq=Δq, Δh=Δh)
end

#####
##### Fixed-point iteration for roughness length
#####

# Replaces `surface_conditions`
@inline function compute_similarity_theory_fluxes(similarity_profile,
                                                  surface_state,
                                                  atmos_state,
                                                  params,
                                                  maxiter=100)
        
    # TODO: Unpack params for 𝜅, 𝑔, thermo_params, h_atmos_boundary_layer, LH_v0
    𝑔 = SFP.grav(params)
    𝜅 = SFP.von_karman_const(params)
    h_atmos_boundary_layer = atmos_state.h_boundary_layer
    LH_v0 = SFP.LH_v0(params)
    thermo_params = params.thermo_params
    FT = eltype(thermo_params)
    turbulent_prandtl_number = 1//3 #FIXME: Get from ClimaParams
    # Initial guess for the similarity scales u★, θ★, q★.
    # Does not really matter if we are sophisticated or not, it converges 
    # in about 10 iterations no matter what...
    u★ = FT(1e-4)
    Σ★ = SimilarityScales(u★, u★, u★) 

    (; Δu, Δv, Δθ, Δq, Δh) = state_differences(surface_state, atmos_state, Σ★, params)

    # The inital velocity scale assumes that
    # the gustiness velocity `Uᴳ` is equal to 0.5 ms⁻¹. 
    # That will be refined later on.
    Uᴳᵢ² = convert(FT, atmos_state.gustiness_parameter^2)
    ΔU = sqrt(Δu^2 + Δv^2 + Uᴳᵢ²)

    # Initialize the solver
    iteration = 0
    Σ₀ = Σ★
    while iterating(Σ★ - Σ₀, iteration, maxiter, similarity_profile)
        Σ₀ = Σ★
        # Refine both the similarity scale and the effective
        # velocity difference ΔU, including gustiness.
        Σ★, ΔU = refine_similarity_variables(Σ★, 
                                             ΔU, 
                                             similarity_profile,
                                             surface_state,
                                             atmos_state,
                                             params) 
        iteration += 1
    end

    u★ = Σ★.momentum
    θ★ = Σ★.temperature
    q★ = Σ★.water_vapor

    # Compute updated state differences
    (; Δu, Δv, Δθ, Δq, Δh) = state_differences(surface_state, atmos_state, Σ★, params)

    θ★ /= turbulent_prandtl_number # TODO: get from ClimaParams
    q★ /= turbulent_prandtl_number # TODO: get from ClimaParams
    
    # `u★² ≡ sqrt(τx² + τy²)`
    # We remove the gustiness by dividing by `ΔU`
    τx = - u★^2 * Δu / ΔU
    τy = - u★^2 * Δv / ΔU

    atmos_ρ = atmos_state.args.ρ
    atmos_ts = TD.PhaseEquil_ρθq(thermo_params,atmos_ρ, atmos_state.θ_a, atmos_state.q_a)
    cp_m = TD.cp_m(thermo_params, atmos_ts) # moist heat capacity

    # Estimate ρ_s using adiabatic extrapolation
    ρ_a = TD.air_density(thermo_params, atmos_ts)
    cv_m = TD.cv_m(thermo_params, TD.PhasePartition(atmos_state.q_a))
    R_m =  TD.gas_constant_air(thermo_params, atmos_ts)
    θ_a = atmos_state.θ_a
    # TODO Define `surface_args`
    #θ_s = surface_variable(surface_state.θ_s, 
    #                       surface_args, 
    #                       similarity_scales, 
    #                       atmos_state, 
    #                       params)
    θ_s = surface_state.θ_s

    ρ_s = ρ_a * (θ_s/θ_a)^(cv_m/R_m)

    # These currently are not consistent with the dycore paper
    # Use ρ_s, use ΔDSE, use LH_v0
    fluxes = (;
              sensible_heat = - ρ_s * cp_m * u★ * θ★, # This needs to be the same as -ρ_s Cd |u|_a ΔDSE
              latent_heat   = - ρ_s * u★ * q★ * LH_v0, # This might be OK
              water_vapor   = - ρ_s * u★ * q★, # This might be OK
              x_momentum    = + ρ_s * τx,
              y_momentum    = + ρ_s * τy,
              r_ae = Δq/(ρ_a * u★ * q★), # Land needs this, and it is not computable internally to land from only the fluxes in coupled simulation.
              scale_vars = (u★, θ★, q★),
    )

    return fluxes
end

@inline function iterating(Σ★, iteration, maxiter, similarity_profile)
    hasnt_started = iteration == 0
    tolerance = sqrt(eps(FT))
    converged = norm(Σ★) < tolerance
    reached_maxiter = iteration ≥ maxiter
    return !(converged | reached_maxiter) | hasnt_started
end

"""
    TODO: Refactor OCEAN model expressions to match `SurfaceStates` 
"""
@inline function refine_similarity_variables(Σ_est, 
                                             ΔU_est,
                                             similarity_profile,
                                             surface_state,
                                             atmos_state,
                                             params)

    𝑔 = SFP.grav(params)
    𝜅 = SFP.von_karman_const(params)
    h_atmos_boundary_layer = atmos_state.h_boundary_layer
    LH_v0 = SFP.LH_v0(params)
    thermo_params = params.thermo_params
    FT = eltype(thermo_params)
    turbulent_prandtl_number = 1//3 #FIXME: Get from ClimaParams
    gustiness = atmos_state.gustiness_parameter
    
    # Update the state differences given the new guess for Σ
    (; Δu, Δv, Δθ, Δq, Δh) = state_differences(surface_state, atmos_state, Σ_est, params)

    # Unpack and compute roughness lengths according to the surface model functions
    (; 𝑧0m, 𝑧0θ, 𝑧0q) = surface_state.roughness_lengths

    # ??
    #z₀q = surface_variable(roughness_length_q,surface_args,similarity_scales,atmos_state,params)
    #z₀b = surface_variable(roughness_length_θ,surface_args,similarity_scales,atmos_state,params)
    #z₀u = surface_variable(roughness_length_mom,surface_args,similarity_scales,atmos_state,params)

    # "initial" scales because we will recompute them
    u★ = Σ_est.momentum
    θ★ = Σ_est.temperature
    q★ = Σ_est.water_vapor
    uτ = ΔU_est
    atmos_ρ = atmos_state.args.ρ
    atmos_ts = TD.PhaseEquil_ρθq(thermo_params,atmos_ρ, atmos_state.θ_a, atmos_state.q_a)

    # Compute Monin-Obukhov length scale depending on a `buoyancy flux`
    b★ = buoyancy_scale(θ★, q★, atmos_ts, thermo_params, 𝑔)
    # Monin-Obhukov similarity length scale and non-dimensional height
    L★ = ifelse(b★ == 0, zero(b★), - u★^3 * atmos_state.θ_a / (u★ * θ★ * 𝜅 * 𝑔))
    ζ = Δh / L★ 

    ψm = UF.psi(ufunc, ζ, UF.MomentumTransport())
    ψs = UF.psi(ufunc, ζ, UF.HeatTransport()) # TODO Rename HeatTransport > ScalarTransport
    ψm₀ = UF.psi(ufunc, 𝑧0m * ζ / Δstate.Δh, UF.MomentumTransport())
    ψh₀ = UF.psi(ufunc, 𝑧0θ * ζ / Δstate.Δh, UF.HeatTransport())
    ψq₀ = UF.psi(ufunc, 𝑧0q(u★,ζ) * ζ / Δstate.Δh, UF.HeatTransport())
 
    # compute rhs in Δχ/u★ = (f(ζ,𝑧0...))
    F_m = log(Δstate.Δh / 𝑧0m) - ψm + ψm₀
    F_h = log(Δstate.Δh / 𝑧0θ) - ψs + ψh₀
    F_q = log(Δstate.Δh / 𝑧0q(u★, ζ)) - ψs + ψq₀

    # Review against nishizawa notation
    χu = 𝜅/F_m 
    χθ = 𝜅/F_m
    χq = 𝜅/F_m

    u★ = χu * uτ
    θ★ = χθ * Δstate.Δθ
    q★ = χq * Δstate.Δq
    
    # Buoyancy flux similarity scale for gustiness (Edson 2013)
    hᵢ = h_atmos_boundary_layer
    Jᵇ = - u★ * b★
    Uᴳ = gustiness * cbrt(Jᵇ * hᵢ)

    # New velocity difference accounting for gustiness
    ΔU = sqrt(Δu^2 + Δv^2 + Uᴳ^2)
    return SimilarityScales(u★, θ★, q★), ΔU
end
