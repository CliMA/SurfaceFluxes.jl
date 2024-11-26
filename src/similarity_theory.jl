### REVAMP 

struct SimilarityScales{U, T, Q}
    momentum :: U
    temperature :: T
    water_vapor :: Q
end

function -(a::SimilarityScales, b::SimilarityScales)
    Δu = a.momentum - b.momentum
    Δθ = a.temperature - b.temperature
    Δq = a.water_vapor - b.water_vapor
    return SimilarityScales(Δu, Δθ, Δq)
end

struct SurfaceState{RL, 
                    U <: Union{Function, FT}, V <: Union{Function, FT}, 
                    Q <: Union{Function, FT}, T <: Union{Function, FT}, 
                    X, FT} <: AbstractStateVariables{FT <: AbstractFloat} 
    roughness_lengths::RL # Should be a named tuple with vapor, momentum, and buoyancy. Those should then be of type Union{Function, FT}
    u_s::U
    v_s::V
    h::FT
    q_s::Q
    θ_s::T
    args::X
    gustiness::FT
end

"""
    state_differences(surface_state, atmos_state, similarity_scales, params)
returns `NamedTuple` of state differences across first interior cell
"""
function state_differences(surface_state, atmos_state, similarity_scales, params)
    q_a = atmos_state.q_a
    surface_args = surface_state.args
    q_s = surface_state.q_s(surface_args, similarity_scales, atmos_state, params)
    Δq = q_a - q_s
    
    u_a = atmos_state.u_a
    u_s = surface_state.u_s(surface_args, similarity_scales, atmos_state, params)
    Δu = u_a - u_s
    
    v_a = atmos_state.v_a
    v_s = surface_state.v_s(surface_args, similarity_scales, atmos_state, params)
    Δv = v_a - v_s
    
    θ_a = atmos_state.θ_a
    θ_s = surface_state.θ_s(surface_args, similarity_scales, atmos_state, params)
    Δθ= θ_a - θ_s
    
    h_a = atmos_state.h_a
    h_s = surface_state.h_s
    Δh = h_a - h_s 

    return (; Δu=Δu, Δv=Δv, Δθ=Δθ, Δq=Δq, Δh=Δh)
end

Statistics.norm(a::SimilarityScales) = norm(a.momentum) + norm(a.temperature) + norm(a.water_vapor)

#####
##### Fixed-point iteration for roughness length
#####

@inline function compute_similarity_theory_fluxes(similarity_theory,
                                                  surface_state,
                                                  atmos_state,
                                                  params,
                                                  maxiter)
        
    # TODO: Unpack params for 𝜅, grav, thermo_params, h_atmos_boundary_layer

    # Prescribed difference between two states
    ℂₐ = thermodynamics_parameters

    # Initial guess for the similarity scales u★, θ★, q★.
    # Does not really matter if we are sophisticated or not, it converges 
    # in about 10 iterations no matter what...
    u★ = convert(eltype(Δh), 1e-4)
    Σ★ = SimilarityScales(u★, u★, u★) 

    Δstate = state_differences(surface_state, atmos_state, Σ★, params)

    # The inital velocity scale assumes that
    # the gustiness velocity `Uᴳ` is equal to 0.5 ms⁻¹. 
    # That will be refined later on.
    FT = eltype(Δstate.Δh)
    Uᴳᵢ² = convert(FT, 0.5^2)
    ΔU = sqrt(Δstate.Δu^2 + Δstate.Δv^2 + Uᴳᵢ²)

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

    θ★ = θ★ / similarity_theory.turbulent_prandtl_number
    q★ = q★ / similarity_theory.turbulent_prandtl_number

    # `u★² ≡ sqrt(τx² + τy²)`
    # We remove the gustiness by dividing by `ΔU`
    τx = - u★^2 * Δu / ΔU
    τy = - u★^2 * Δv / ΔU

    𝒬ₐ = atmos_state.ts
    ρₐ = AtmosphericThermodynamics.air_density(ℂₐ, 𝒬ₐ)
    cₚ = AtmosphericThermodynamics.cp_m(ℂₐ, 𝒬ₐ) # moist heat capacity
    ℰv = AtmosphericThermodynamics.latent_heat_vapor(ℂₐ, 𝒬ₐ)

    fluxes = (;
        sensible_heat = - ρₐ * cₚ * u★ * θ★,
        latent_heat   = - ρₐ * u★ * q★ * ℰv,
        water_vapor   = - ρₐ * u★ * q★,
        x_momentum    = + ρₐ * τx,
        y_momentum    = + ρₐ * τy,
    )

    return fluxes
end

"""
    TODO: Refactor OCEAN model expressions to match `SurfaceStates` 
"""
@inline function refine_similarity_variables(estimated_similarity_scales, 
                                             velocity_scale,
                                             similarity_theory,
                                             surface_state,
                                             atmos_state,
                                             params)
    # TODO Check param unpack
    
    Δstate = state_differences(surface_state, atmos_state, estimated_similarity_scales, params)
    h  = Δstate.Δh
    (z0m, z0θ, z0b) = surface_state.roughness_lengths
    
    ζ = # Given u★ = u★²/𝜅b★ # Check signs
    
    z₀b = surface_state.roughness_lengths.z₀b(surface_args, similarity_scales, atmos_state, params)
    z₀u = surface_state.roughness_lengths.z₀u(surface_args, similarity_scales, atmos_state, params)

    # Current SurfaceFluxes calls to stability functions 
    ψₘ = UF.psi(similarity_theory, ζ, UF.MomentumTransport())
    ψₕ = UF.psi(similarity_theory, ζ, UF.ScalarTransport())
    ψₘ₀ = UF.psi(similarity_theory, z₀u * ζ / h, UF.MomentumTransport())
    ψₕ₀ = UF.psi(similarity_theory, z₀b * ζ / h, UF.ScalarTransport())
    
    denominator = log(Δh / z0(sc, transport)) - ψ + ψ₀

    # "initial" scales because we will recompute them
    u★ = estimated_similarity_scales.momentum
    θ★ = estimated_similarity_scales.temperature
    q★ = estimated_similarity_scales.water_vapor
    uτ = velocity_scale

    # Extract roughness lengths
    gustiness = surface_states.gustiness_parameter 
    thermo_params  = thermodynamics_parameters
    g  = gravitational_acceleration

    # Compute Monin-Obukhov length scale depending on a `buoyancy flux`
    # TODO: Check buoyancy_scale function definitions and conventions
    b★ = buoyancy_scale(θ★, q★, thermo_params)
    # Monin-Obhukov similarity length scale and non-dimensional height
    ϰ  = von_karman_constant
    L★ = ifelse(b★ == 0, zero(b★), - u★^2 / (ϰ * b★))
    
    # Transfer coefficients at height `h`
    profile_type = similarity_theory.similarity_profile_type
    χu = ϰ /denom
    χθ = ϰ /denom
    χq = ϰ /denom

    Δu = differences.u
    Δv = differences.v
    Δθ = differences.θ
    Δq = differences.q

    # u★ including gustiness
    u★ = χu * uτ
    θ★ = χθ * Δθ
    q★ = χq * Δq

    # Buoyancy flux similarity scale for gustiness (Edson 2013)
    hᵢ = atmos_boundary_layer_height
    Jᵇ = - u★ * b★
    Uᴳ = gustiness * cbrt(Jᵇ * hᵢ)

    # New velocity difference accounting for gustiness
    ΔU = sqrt(Δu^2 + Δv^2 + Uᴳ^2)

    return SimilarityScales(u★, θ★, q★), ΔU
end
