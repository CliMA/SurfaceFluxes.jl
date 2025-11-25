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

import Thermodynamics
const TD = Thermodynamics

import .UniversalFunctions
const UF = UniversalFunctions

import .Parameters

const SFP = Parameters
const APS = SFP.AbstractSurfaceFluxesParameters

include("types.jl")
include("utilities.jl")
include("physical_scale_coefficient_methods.jl")
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

The main user-facing function of the module. Computes surface conditions
based on Monin-Obukhov similarity theory.

## Arguments
- `param_set`: Parameter set containing physical and thermodynamic constants
- `sc`: Surface conditions container (Fluxes, ValuesOnly, Coefficients, or FluxesAndFrictionVelocity)
- `scheme`: Discretization scheme (PointValueScheme for finite difference or LayerAverageScheme for finite volume)
- `tol_neutral`: Tolerance for neutral stability detection based on `ΔDSEᵥ` (default: `cp_d / 100`)
- `tol`: Convergence tolerance for iterative solver (default: `sqrt(eps(FT))`)
- `maxiter`: Maximum number of iterations (default: 10)

## Returns
Returns a `SurfaceFluxConditions` struct containing:
  - `L_MO`:   Monin-Obukhov lengthscale [m]
  - `shf`:    Sensible heat flux [W/m²]
  - `lhf`:    Latent heat flux [W/m²]
  - `ρτxz`:   Momentum flux, eastward component [kg/(m·s²)]
  - `ρτyz`:   Momentum flux, northward component [kg/(m·s²)]
  - `ustar`:  Friction velocity [m/s]
  - `Cd`:     Momentum exchange coefficient
  - `Ch`:     Heat exchange coefficient
  - `E`:      Evaporation rate [kg/(m²·s)]
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
    tol::FT = eps(FT),
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

    obukhov_similarity_solution(
        param_set::AbstractSurfaceFluxesParameters,
        sc::AbstractSurfaceConditions,
        scheme,
        tol,
        tol_neutral,
        maxiter,
    )

Compute and return the Monin-Obukhov similarity solution.

Solves for the Monin-Obukhov lengthscale (L_MO) and related similarity scales
using an iterative Newton-Raphson method. The solution depends on the
particular surface condition type `sc <: AbstractSurfaceConditions`.

## Arguments
- `param_set`: Parameter set containing physical constants
- `sc`: Surface conditions container
- `scheme`: Discretization scheme
- `tol`: Convergence tolerance for iterative solver
- `tol_neutral`: Tolerance for neutral stability detection
- `maxiter`: Maximum number of iterations

## Returns
Returns a named tuple containing:
  - `L★`: Monin-Obukhov lengthscale [m]
  - `u★`: Friction velocity [m/s]
  - `DSEᵥ★`: Virtual dry static energy scale [J/kg]
  - `θᵥ★`: Virtual potential temperature scale [K]
  - `q★`: Specific humidity scale [kg/kg]
  - `𝓁u`: Momentum roughness length [m]
  - `𝓁θ`: Heat roughness length [m]
  - `𝓁q`: Moisture roughness length [m]
"""
function obukhov_similarity_solution end

obukhov_similarity_solution(sfc::SurfaceFluxConditions) = sfc.L_MO

function compute_Fₘₕ(sc, ufₛ, ζ, 𝓁, transport)
    return log(Δz(sc) / 𝓁) -
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
        X★₀ = (u★ = u★₀, DSEᵥ★ = FT(δ), θᵥ★ = FT(δ), q★ = FT(δ),
            L★ = FT(10),
            𝓁u = 𝓁u₀, 𝓁θ = 𝓁θ₀, 𝓁q = 𝓁q₀)
        X★ = obukhov_iteration(X★₀, sc, scheme, param_set, tol, tol_neutral)
        return X★
    else
        X★₀ = (u★ = u★₀, DSEᵥ★ = FT(δ), θᵥ★ = FT(δ), q★ = FT(δ),
            L★ = FT(-10),
            𝓁u = 𝓁u₀, 𝓁θ = 𝓁θ₀, 𝓁q = 𝓁q₀)
        X★ = obukhov_iteration(X★₀, sc, scheme, param_set, tol, tol_neutral)
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
    iterate_interface_fluxes()
"""
function iterate_interface_fluxes(sc::Union{ValuesOnly, Fluxes},
    q_surface,
    approximate_interface_state,
    atmosphere_state,
    surface_state,
    scheme::SolverScheme,
    param_set::APS,
    tol_neutral,
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
    Tₛ = surface_temperature(param_set, sc, (; u★, q★))

    # Surface Quantities and state differences
    surface_args = sc.state_sfc.args
    ΔU = windspeed(sc)

    ### Compute Monin--Obukhov length scale depending on the buoyancy scale b★
    ### The windspeed function accounts for a wind-gust parameter.
    b★ = DSEᵥ★ * 𝑔 / DSEᵥ_in(param_set, sc)
    L★ = ifelse(abs(ΔDSEᵥ(param_set, sc)) <= tol_neutral,
        sign(ΔDSEᵥ(param_set, sc)) * FT(Inf),
        non_zero(u★^2 / (𝜅 * b★)))
    ## The new L★ estimate is then used to update all scale variables
    ## with stability correction functions (compute_Fₘₕ)
    ζ = non_zero(Δz(sc) / L★)

    ### Compute new values for the scale parameters given the relation
    χu = 𝜅 / compute_Fₘₕ(sc, uf, ζ, 𝓁u, UF.MomentumTransport())
    χDSEᵥ = 𝜅 / compute_Fₘₕ(sc, uf, ζ, 𝓁θ, UF.HeatTransport())
    χq = 𝜅 / compute_Fₘₕ(sc, uf, ζ, 𝓁q, UF.HeatTransport())
    χθᵥ = 𝜅 / compute_Fₘₕ(sc, uf, ζ, 𝓁θ, UF.HeatTransport())

    ## Re-compute scale variables
    u★ = χu * ΔU
    DSEᵥ★ = χDSEᵥ * ΔDSEᵥ(param_set, sc)
    q★ = χq * Δq
    θᵥ★ = χθᵥ * Δθᵥ(param_set, sc)

    # Returns a NamedTuple with similarity scales and roughness lengths:
    # - u★: Friction velocity [m/s]
    # - DSEᵥ★: Virtual dry static energy scale [J/kg]
    # - q★: Specific humidity scale [kg/kg]
    # - L★: Monin-Obukhov lengthscale [m]
    # - θᵥ★: Virtual potential temperature scale [K]
    # - 𝓁u: Momentum roughness length [m]
    # - 𝓁θ: Heat roughness length [m]
    # - 𝓁q: Moisture roughness length [m]
    return (; u★, DSEᵥ★, q★, L★, θᵥ★, 𝓁u, 𝓁θ, 𝓁q)
end

function obukhov_iteration(X★,
    sc,
    scheme,
    param_set,
    tol,
    tol_neutral,
    maxiter = 40,
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
            param_set,
            tol_neutral)
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
