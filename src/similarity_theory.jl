### REVAMP 

"""
    SimilarityScales

Holds the Monin-Obukhov solution state u★, θ★, q★
"""
struct SimilarityScales{U, T, Q}
    momentum::U # ustar
    temperature::T # bstar or theta-star?
    water_vapor::Q # qstar
end

"""
    -(a::SimilarityScales, b::SimilarityScales)

Component-wise subtraction
"""
function -(a::SimilarityScales, b::SimilarityScales)
    Δu = a.momentum - b.momentum
    Δb = a.temperature - b.temperature
    Δq = a.water_vapor - b.water_vapor
    return SimilarityScales(Δu, Δθ, Δq)
end

Statistics.norm(a::SimilarityScales) = norm(a.momentum) + norm(a.temperature) + norm(a.water_vapor)

"""
    SurfaceStates

Holds the functions required to compute roughness lengths, surface conditions
"""
struct SurfaceState{RL, 
                    U <: Union{Function, FT},
                    V <: Union{Function, FT}, 
                    Q <: Union{Function, FT},
                    T <: Union{Function, FT}, 
                    X,
                    FT} <: AbstractStateVariables where {FT <: AbstractFloat} 
    roughness_lengths::RL # Should be a named tuple with vapor, momentum, and buoyancy. Those should then be of type Union{Function, FT}
    u_s::U
    v_s::V
    q_s::Q
    θ_s::T
    args::X # whatever is needed to compute surface conditions, provided by the surface model.
    h::FT
    gustiness::FT
end

function surface_variable(surface_state_var::Function, surface_args, similarity_scales, atmos_state, params)
    return surface_state_var(surface_args, similarity_scales, atmos_state, params)
end

function surface_variable(surface_state_var::FT, _...) where {FT <: AbstractFloat}
    return surface_state_var
end
    
"""
    state_differences(surface_state, atmos_state, similarity_scales, params)

returns `NamedTuple` of state differences across first interior cell
"""
function state_differences(surface_state, atmos_state, similarity_scales, params)
    surface_args = surface_state.args
    
    q_a = atmos_state.q_a
    q_s = surface_variable(surface_state.q_s, surface_args, similarity_scales, atmos_state, params)
    Δq = q_a - q_s
    
    u_a = atmos_state.u_a
    u_s = surface_variable(surface_state.u_s, surface_args, similarity_scales, atmos_state, params)
    Δu = u_a - u_s
    
    v_a = atmos_state.v_a
    v_s = surface_variable(surface_state.v_s, surface_args, similarity_scales, atmos_state, params)
    Δv = v_a - v_s
    
    θ_a = atmos_state.θ_a
    θ_s = surface_variable(surface_state.θ_s, surface_args, similarity_scales, atmos_state, params)
    Δθ= θ_a - θ_s
    
    h_a = atmos_state.h_a
    h_s = surface_state.h_s
    Δh = h_a - h_s 

    return (; Δu=Δu, Δv=Δv, Δθ=Δθ, Δq=Δq, Δh=Δh)
end



#####
##### Fixed-point iteration for roughness length
#####

@inline function compute_similarity_theory_fluxes(similarity_theory,
                                                  surface_state,
                                                  atmos_state,
                                                  params,
                                                  maxiter)
        
    # TODO: Unpack params for 𝜅, grav, thermo_params, h_atmos_boundary_layer
    (; grav, thermo_params, 𝜅, h_atmos_boundary_layer) = params
    FT = typeof(grav)
    # Initial guess for the similarity scales u★, θ★, q★.
    # Does not really matter if we are sophisticated or not, it converges 
    # in about 10 iterations no matter what...
    u★ = FT(1e-4)
    Σ★ = SimilarityScales(u★, u★, u★) 

    (; Δu, Δv, Δθ, Δq, Δh) = state_differences(surface_state, atmos_state, Σ★, params)

    # The inital velocity scale assumes that
    # the gustiness velocity `Uᴳ` is equal to 0.5 ms⁻¹. 
    # That will be refined later on.
    Uᴳᵢ² = convert(FT, 0.5^2)
    ΔU = sqrt(Δu^2 + Δv^2 + Uᴳᵢ²)

    # Initialize the solver
    iteration = 0
    Σ₀ = Σ★
    while iterating(Σ★ - Σ₀, iteration, maxiter, similarity_theory)
        Σ₀ = Σ★
        # Refine both the similarity scale and the effective
        # velocity difference ΔU, including gustiness.
        Σ★, ΔU = refine_similarity_variables(Σ★, 
                                             ΔU, 
                                             similarity_theory,
                                             surface_state,
                                             atmos_state,
                                             params) 
        # params contains 
        #       thermo_params
        #       atmos boundary layer height
        #       von karman constant
        #       gravitational acceleration
        iteration += 1
    end

    u★ = Σ★.momentum
    θ★ = Σ★.temperature
    q★ = Σ★.water_vapor

    # Compute updated state differences
    (; Δu, Δv, Δθ, Δq, Δh) = state_differences(surface_state, atmos_state, Σ★, params)

    θ★ = θ★ / similarity_theory.turbulent_prandtl_number
    q★ = q★ / similarity_theory.turbulent_prandtl_number
    
    # `u★² ≡ sqrt(τx² + τy²)`
    # We remove the gustiness by dividing by `ΔU`
    τx = - u★^2 * Δu / ΔU
    τy = - u★^2 * Δv / ΔU

    atmos_ts = atmos_state.ts
    ρ_a= AtmosphericThermodynamics.air_density(thermo_params, atmos_ts)
    cp_m = AtmosphericThermodynamics.cp_m(thermo_params, atmos_ts) # moist heat capacity
    LH_v= AtmosphericThermodynamics.latent_heat_vapor(thermo_params, atmos_ts)

    # These currently are not consistent with the dycore paper
    # Use ρ_s, use ΔDSE (use LH_v0?)
    fluxes = (;
              sensible_heat = - ρ_a * cp_m * u★ * θ★,
              latent_heat   = - ρ_a * u★ * q★ * LH_v,
              water_vapor   = - ρ_a * u★ * q★,
              x_momentum    = + ρ_a * τx,
              y_momentum    = + ρ_a * τy,
              r_ae = Δq/(ρ_a * u★ * q★) # Land needs this, and it is not computable internally to land from only the fluxes in coupled simulation.
    )

    return fluxes
end

"""
    TODO: Refactor OCEAN model expressions to match `SurfaceStates` 
"""
@inline function refine_similarity_variables(Σ_est, 
                                             ΔU_est,
                                             similarity_theory,
                                             surface_state,
                                             atmos_state,
                                             params)

    (; grav, thermo_params, 𝜅, h_atmos_boundary_layer) = params
    gustiness = surface_states.gustiness_parameter
    
    # Update the state differences given the new guess for Σ
    (; Δu, Δv, Δθ, Δq, Δh) = state_differences(surface_state, atmos_state, Σ_est, params)

    # Unpack and compute roughness lengths according to the surface model functions
    (; roughness_length_mom, roughness_length_θ, roughness_length_q) = surface_state.roughness_lengths

    z₀q = surface_variable(roughness_length_q, surface_args, similarity_scales, atmos_state, params)
    z₀b = surface_variable(roughness_length_θ, surface_args, similarity_scales, atmos_state, params)
    z₀u = surface_variable(roughness_length_mom, surface_args, similarity_scales, atmos_state, params)

    # "initial" scales because we will recompute them
    u★ = Σ_est.momentum
    θ★ = Σ_est.temperature
    q★ = Σ_est.water_vapor
    uτ = ΔU_est

    # Compute Monin-Obukhov length scale depending on a `buoyancy flux`
    # TODO: Check buoyancy_scale function definitions and conventions
    b★ = buoyancy_scale(θ★, q★, thermo_params)
    # Monin-Obhukov similarity length scale and non-dimensional height
    L★ = ifelse(b★ == 0, zero(b★), - u★^2 / (ϰ * b★))

    # Current SurfaceFluxes calls to stability functions 
    
    ζ = # Given u★, L = u★²/𝜅b★ # Check signs
    ψₘ = UF.psi(similarity_theory, ζ, UF.MomentumTransport())
    ψₕ = UF.psi(similarity_theory, ζ, UF.ScalarTransport())
    ψₘ₀ = UF.psi(similarity_theory, z₀u * ζ / h, UF.MomentumTransport())
    ψₕ₀ = UF.psi(similarity_theory, z₀b * ζ / h, UF.ScalarTransport())
    
    denom = log(Δh / z0(sc, transport)) - ψ + ψ₀

    # Transfer coefficients at height `h`
    χu = ϰ /denom
    χθ = ϰ /denom
    χq = ϰ /denom

    # u★ including gustiness
    u★ = χu * uτ
    θ★ = χθ * Δθ
    q★ = χq * Δq

    # Buoyancy flux similarity scale for gustiness (Edson 2013)
    hᵢ = h_atmos_boundary_layer
    Jᵇ = - u★ * b★
    Uᴳ = gustiness * cbrt(Jᵇ * hᵢ)

    # New velocity difference accounting for gustiness
    ΔU = sqrt(Δu^2 + Δv^2 + Uᴳ^2)

    return SimilarityScales(u★, θ★, q★), ΔU
end
