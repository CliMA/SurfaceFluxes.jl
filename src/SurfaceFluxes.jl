"""
    SurfaceFluxes

## Interface
  - [`surface_conditions`](@ref) computes
    - Monin-Obukhov length
    - Potential temperature flux (if not given) using Monin-Obukhov theory
    - transport fluxes using Monin-Obukhov theory
    - friction velocity/temperature scale/tracer scales
    - exchange coefficients

## References
 - [Nishizawa2018](@cite)
 - [Byun1990](@cite)

"""
module SurfaceFluxes
import RootSolvers
const RS = RootSolvers
include("UniversalFunctions.jl")
include("Parameters.jl")

using DocStringExtensions
const DSE = DocStringExtensions

import Thermodynamics
const TD = Thermodynamics

import .UniversalFunctions
const UF = UniversalFunctions

import .Parameters

const SFP = Parameters
const APS = SFP.AbstractSurfaceFluxesParameters

include("types.jl")
include("utilities.jl")
include("roughness_models.jl")
include("coefficient_inputs.jl")
include("evaporation_methods.jl")
include("latent_heat_methods.jl")
include("sensible_heat_methods.jl")
include("momentum_exchange_coefficient_methods.jl")
include("heat_exchange_coefficient_methods.jl")
include("friction_velocity_methods.jl")
include("buoyancy_flux_methods.jl")
include("profile_recovery.jl")

"""
    surface_conditions(
        param_set::AbstractSurfaceFluxesParameters,
        sc::SurfaceFluxes.AbstractSurfaceConditions,
        scheme::SurfaceFluxes.SolverScheme = PointValueScheme();
        tol_neutral = SFP.cp_d(param_set) / 100,
        tol = sqrt(eps(FT)),
        maxiter::Int = 10,
    )

The main user facing function of the module.
It computes the surface conditions
based on the Monin-Obukhov similarity functions. Requires
information about thermodynamic parameters (`param_set`),
the surface state `sc`, and the discretisation `scheme`. Default tolerance for
Monin-Obukhov length is absolute (i.e. has units [m]).
Returns the RootSolvers `CompactSolution` by default.

Result struct of type SurfaceFluxConditions contains:
  - L_MO:   Monin-Obukhov lengthscale
  - shf:    Sensible Heat Flux
  - lhf:    Latent Heat Flux
  - ρτxz:   Momentum Flux (Eastward component)
  - ρτyz:   Momentum Flux (Northward component)
  - ustar:  Friction velocity
  - Cd:     Momentum Exchange Coefficient
  - Ch:     Thermal Exchange Coefficient
  - z₀:     Aerodynamic roughness lengths
"""
function surface_conditions(
    param_set::APS{FT},
    sc::AbstractSurfaceConditions,
    scheme::SolverScheme = PointValueScheme();
    tol_neutral = SFP.cp_d(param_set) / 100,
    tol = sqrt(eps(FT)),
    maxiter::Int = 30,
) where {FT}
    uft = SFP.uf_params(param_set)
    X★ = obukhov_similarity_solution(
        param_set,
        sc,
        scheme,
        tol,
        tol_neutral,
        maxiter,
    )
    L_MO = X★.L★
    ustar = X★.u★
    𝓁u = X★.𝓁u
    𝓁θ = X★.𝓁θ
    Cd = momentum_exchange_coefficient(
        param_set,
        L_MO,
        ustar,
        sc,
        scheme,
        tol_neutral,
    )
    Ch =
        heat_exchange_coefficient(
            param_set,
            L_MO,
            ustar,
            sc,
            scheme,
            tol_neutral,
        )
    shf = sensible_heat_flux(param_set, Ch, sc, scheme)
    lhf = latent_heat_flux(param_set, Ch, sc, scheme)
    buoy_flux = compute_buoyancy_flux(
        param_set,
        shf,
        lhf,
        ts_in(sc),
        ts_sfc(sc),
        scheme,
    )
    ρτxz, ρτyz = momentum_fluxes(param_set, Cd, sc, scheme)
    E = evaporation(param_set, sc, Ch)
    return SurfaceFluxConditions(
        L_MO,
        shf,
        lhf,
        buoy_flux,
        ρτxz,
        ρτyz,
        ustar,
        Cd,
        Ch,
        E,
    )
end

function surface_conditions(
    param_set::APS{FT},
    sc::FluxesAndFrictionVelocity,
    scheme::SolverScheme = PointValueScheme();
    tol_neutral = SFP.cp_d(param_set) / 100,
    tol::FT = sqrt(eps(FT)),
    maxiter::Int = 10,
) where {FT}
    uft = SFP.uf_params(param_set)
    X★ = obukhov_similarity_solution(
        param_set,
        sc,
        scheme,
        tol,
        tol_neutral,
        maxiter,
    )
    Cd = momentum_exchange_coefficient(param_set, X★.L★, X★.u★, sc, scheme, tol_neutral)
    Ch = heat_exchange_coefficient(param_set, X★.L★, X★.u★, sc, scheme, tol_neutral)
    shf = sc.shf
    lhf = sc.lhf
    buoy_flux = compute_buoyancy_flux(
        param_set,
        shf,
        lhf,
        ts_in(sc),
        ts_sfc(sc),
        scheme,
    )
    ρτxz, ρτyz = momentum_fluxes(param_set, Cd, sc, scheme)
    E = evaporation(param_set, sc, Ch)
    return SurfaceFluxConditions(
        X★.L★,
        shf,
        lhf,
        buoy_flux,
        ρτxz,
        ρτyz,
        X★.u★,
        Cd,
        Ch,
        E,
    )
end

"""
    obukhov_similarity_solution(sfc::SurfaceFluxConditions)

    obukhov_similarity_solution( # internal method
        param_set::AbstractSurfaceFluxesParameters,
        sc::AbstractSurfaceConditions,
        uft,
        scheme,
        tol,
        tol_neutral,
        maxiter,
    )

Compute and return the Monin-Obukhov lengthscale (LMO).

The internal method for computing LMO depends on the
particular surface condition `sc <: AbstractSurfaceConditions`. 
"""
function obukhov_similarity_solution end

obukhov_similarity_solution(sfc::SurfaceFluxConditions) = sfc.L_MO

function compute_Fₘₕ(sc, ufₛ, ζ, 𝓁, transport)
    return log(Δz(sc)/𝓁) -
           UF.psi(ufₛ, ζ, transport) +
           UF.psi(ufₛ, 𝓁 * ζ / Δz(sc), transport)
end

function obukhov_similarity_solution(
    param_set::APS{FT},
    sc::Union{Fluxes, ValuesOnly},
    scheme,
    tol,
    tol_neutral,
    maxiter,
) where {FT}
    thermo_params = SFP.thermodynamics_params(param_set)
    ufparams = SFP.uf_params(param_set)
    grav = SFP.grav(param_set)
    δ = sign(ΔDSEᵥ(param_set, sc))
    u★₀ = FT(0.1)
    𝓁u₀ = compute_z0(u★₀, param_set, sc, sc.roughness_model, UF.MomentumTransport())
    𝓁θ₀ = compute_z0(u★₀, param_set, sc, sc.roughness_model, UF.HeatTransport())
    𝓁q₀ = compute_z0(u★₀, param_set, sc, sc.roughness_model, UF.HeatTransport())
    # Initial guesses for MOST iterative solution
    if ΔDSEᵥ(param_set, sc) >= FT(0)
        X★₀ = (u★ = u★₀, DSEᵥ★ = FT(δ), θᵥ★=FT(δ), q★ = FT(δ),
            L★ = FT(10),
            𝓁u = 𝓁u₀, 𝓁θ = 𝓁θ₀, 𝓁q = 𝓁q₀)
        X★ = obukhov_iteration(X★₀, sc, scheme, param_set, tol)
        return X★
    else
        X★₀ = (u★ = u★₀, DSEᵥ★ = FT(δ), θᵥ★=FT(δ), q★ = FT(δ),
            L★ = FT(-10),
            𝓁u = 𝓁u₀, 𝓁θ = 𝓁θ₀, 𝓁q = 𝓁q₀)
        X★ = obukhov_iteration(X★₀, sc, scheme, param_set, tol)
        return X★
    end
end

function obukhov_similarity_solution(
    param_set,
    sc::FluxesAndFrictionVelocity,
    scheme,
    args...,
)
    return (L★ = -sc.ustar^3 / SFP.von_karman_const(param_set) /
                 non_zero(compute_buoyancy_flux(param_set, sc, scheme)), u★ = sc.ustar)
end


"""
    momentum_fluxes(param_set, Cd, sc, scheme)

Compute and return the momentum fluxes
## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - Cd: Momentum exchange coefficient
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - scheme: Discretization scheme (currently supports FD and FV)
"""
function momentum_fluxes(param_set, Cd, sc::AbstractSurfaceConditions, scheme)
    thermo_params = SFP.thermodynamics_params(param_set)
    ρ_sfc = TD.air_density(thermo_params, ts_sfc(sc))
    ρτxz = -ρ_sfc * Cd * Δu1(sc) * windspeed(sc)
    ρτyz = -ρ_sfc * Cd * Δu2(sc) * windspeed(sc)
    return (ρτxz, ρτyz)
end

"""
    compute_physical_scale_coeff(param_set, sc, L_MO, transport, ::LayerAverageScheme)

Computes the coefficient for the physical scale of a variable based on Nishizawa(2018)
for the FV scheme.

## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - L_MO: Monin-Obukhov length
  - transport: Transport type, (e.g. Momentum or Heat, used to determine physical scale coefficients)
  - scheme: Discretization scheme (currently supports FD and FV)
"""
function compute_physical_scale_coeff(
    param_set::APS,
    sc::Union{ValuesOnly, Fluxes, FluxesAndFrictionVelocity},
    L_MO,
    𝓁,
    transport,
    ::LayerAverageScheme,
)
    𝜅 = SFP.von_karman_const(param_set)
    uf = SFP.uf_params(param_set)
    π_group = UF.π_group(uf, transport)
    R_z0 = 1 - 𝓁 / Δz(sc)
    denom1 = log(Δz(sc) / 𝓁)
    denom2 = -UF.Psi(uf, Δz(sc) / L_MO, transport)
    denom3 =
        𝓁 / Δz(sc) *
        UF.Psi(uf, 𝓁 / L_MO, transport)
    denom4 = R_z0 * (UF.psi(uf, 𝓁 / L_MO, transport) - 1)
    Σterms = denom1 + denom2 + denom3 + denom4
    return 𝜅 / (π_group * Σterms)
end

"""
    compute_physical_scale_coeff(param_set, sc, L_MO, transport, ::PointValueScheme)

Computes the coefficient for the physical scale of a variable based on Byun (1990)
for the Finite Differences scheme.

## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - L_MO: Monin-Obukhov length
  - transport: Transport type, (e.g. Momentum or Heat, used to determine physical scale coefficients)
  - scheme: Discretization scheme (currently supports FD and FV)
"""
function compute_physical_scale_coeff(
    param_set,
    sc::Union{ValuesOnly, Fluxes, FluxesAndFrictionVelocity},
    L_MO,
    𝓁,
    transport,
    ::PointValueScheme,
)
    𝜅 = SFP.von_karman_const(param_set)
    FT = eltype(𝜅)
    uf = SFP.uf_params(param_set)
    π_group = UF.π_group(uf, transport)
    denom1 = log(FT(Δz(sc) / 𝓁))
    denom2 = -UF.psi(uf, FT(Δz(sc) / L_MO), transport)
    denom3 = UF.psi(uf, FT(𝓁 / L_MO), transport)
    Σterms = denom1 + denom2 + denom3
    return 𝜅 / (π_group * Σterms)
end

@inline function buoyancy_scale(θᵥ★, q★, thermo_params, ts, 𝑔)
    FT = eltype(𝑔)
    Tᵥ = TD.virtual_temperature(thermo_params, ts)
    qₐ = TD.vapor_specific_humidity(thermo_params, ts)
    ε = TD.Parameters.Rv_over_Rd(thermo_params)
    δ = ε - FT(1)
    b★ = 𝑔 / Tᵥ * (θᵥ★ * (1 + δ * qₐ) + δ * Tᵥ * q★)
    return FT(b★)
end

"""
    iterate_interface_fluxes()
"""
function iterate_interface_fluxes(sc::Union{ValuesOnly, Fluxes},
    q_surface,
    approximate_interface_state,
    atmosphere_state,
    surface_state,
    scheme::SolverScheme,
    param_set::APS,
)
    ### Parameter sets
    uf = SFP.uf_params(param_set)
    thermo_params = SFP.thermodynamics_params(param_set)
    𝜅 = SFP.von_karman_const(param_set)
    𝑔 = SFP.grav(param_set)
    FT = eltype(𝑔)
    
    ## "Initial" approximate scales because we will recompute them
    ## Updated values of these will populate the resulting named-tuple
    u★ = approximate_interface_state.u★
    qₛ = qt_sfc(param_set, sc)
    Δq = Δqt(param_set, sc)
    DSEᵥ★ = approximate_interface_state.DSEᵥ★
    θᵥ★ = approximate_interface_state.θᵥ★
    q★ = Δq == eltype(𝑔)(0) ? approximate_interface_state.q★ : eltype(𝑔)(0)
    L★ = approximate_interface_state.L★
    𝓁u = compute_z0(u★, param_set, sc, sc.roughness_model, UF.MomentumTransport())
    𝓁θ = compute_z0(u★, param_set, sc, sc.roughness_model, UF.HeatTransport())
    𝓁q = compute_z0(u★, param_set, sc, sc.roughness_model, UF.HeatTransport())
    Tₛ = surface_temperature(param_set, sc, (;u★, q★))

    # Surface Quantities and state differences
    surface_args = sc.state_sfc.args
    Δdseᵥ = ΔDSEᵥ(param_set, sc)
    ΔU = sqrt(windspeed(sc)^2)

    ### Compute Monin--Obukhov length scale depending on the buoyancy scale b★
    ### The windspeed function accounts for a wind-gust parameter.
    b★ = buoyancy_scale(θᵥ★, q★, thermo_params, ts_sfc(sc), 𝑔)
    L★ = ifelse(b★ == 0, sign(ΔDSEᵥ(param_set, sc)) * FT(Inf), u★^2 / (𝜅 * b★))
    ## The new L★ estimate is then used to update all scale variables
    ## with stability correction functions (compute_Fₘₕ)
    ζ = Δz(sc) / L★

    ### Compute new values for the scale parameters given the relation
    ### Following MOST, χ/χ★ = ψ(ζ, 𝓁, z)
    χu = 𝜅 / compute_Fₘₕ(sc, uf, ζ, 𝓁u, UF.MomentumTransport())
    χDSEᵥ = 𝜅 / compute_Fₘₕ(sc, uf, ζ, 𝓁θ, UF.HeatTransport())
    χq = 𝜅 / compute_Fₘₕ(sc, uf, ζ, 𝓁q, UF.HeatTransport())
    χθᵥ = 𝜅 / compute_Fₘₕ(sc, uf, ζ, 𝓁θ, UF.HeatTransport())

    ## Re-compute scale variables
    u★ = χu * ΔU
    DSEᵥ★ = χDSEᵥ * ΔDSEᵥ(param_set, sc) 
    q★ = χq * Δq
    θᵥ★ = χθᵥ * Δθᵥ(param_set, sc)

    return (;u★, DSEᵥ★, q★, L★, θᵥ★, 𝓁u, 𝓁θ, 𝓁q)
end

function obukhov_iteration(X★,
    sc,
    scheme,
    param_set,
    tol,
    maxiter = 10,
)
    FT = eltype(X★)
    qₛ = surface_specific_humidity(param_set, sc)
    for iter in 1:maxiter
        X★₀ = X★
        X★ = iterate_interface_fluxes(sc,
            qₛ,
            X★₀,
            ts_in(sc),
            ts_sfc(sc),
            scheme,
            param_set)
        if abs(X★.L★ - X★₀.L★) ≤ tol &&
           abs(X★.u★ - X★₀.u★) ≤ tol &&
           abs(X★.q★ - X★₀.q★) ≤ tol &&
           abs(X★.DSEᵥ★ - X★₀.DSEᵥ★) ≤ tol
            break
        end
    end
    return X★
end

# For backwards compatibility with package extensions
if !isdefined(Base, :get_extension)
    include(joinpath("..", "ext", "CreateParametersExt.jl"))
end

end # SurfaceFluxes module
