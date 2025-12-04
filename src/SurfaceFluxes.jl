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
- `tol_neutral`: Tolerance for neutral stability detection (unused, kept for backward compatibility)
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
    # Default to FixedPointIteration, but can be overridden
    solver_method = FixedPointIteration()
    if ΔDSEᵥ(param_set, sc) >= FT(0)
        X★₀ = (u★ = u★₀, DSEᵥ★ = FT(δ), θᵥ★ = FT(δ), q★ = FT(δ),
            L★ = FT(10),
            𝓁u = 𝓁u₀, 𝓁θ = 𝓁θ₀, 𝓁q = 𝓁q₀)
        X★ = obukhov_iteration(X★₀, sc, scheme, param_set, tol, maxiter, solver_method)
        return X★
    else
        X★₀ = (u★ = u★₀, DSEᵥ★ = FT(δ), θᵥ★ = FT(δ), q★ = FT(δ),
            L★ = FT(-10),
            𝓁u = 𝓁u₀, 𝓁θ = 𝓁θ₀, 𝓁q = 𝓁q₀)
        X★ = obukhov_iteration(X★₀, sc, scheme, param_set, tol, maxiter, solver_method)
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
    L★ = u★^2 / (𝜅 * b★)
    ## The new L★ estimate is then used to update all scale variables
    ## with stability correction functions (compute_Fₘₕ)
    ζ = Δz(sc) / L★

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

"""
    bulk_richardson_number_target(param_set, sc)

Compute the target bulk Richardson number from input state variables.

The bulk Richardson number is defined as:
    Ri_b = (g z (DSE_v,in - DSE_v,sfc)) / (DSE_v,sfc |u(z)|^2)
"""
function bulk_richardson_number_target(param_set, sc)
    𝑔 = SFP.grav(param_set)
    FT = eltype(Δz(sc))
    z = z_in(sc)
    DSEᵥ_in_val = DSEᵥ_in(param_set, sc)
    DSEᵥ_sfc_val = DSEᵥ_sfc(param_set, sc)
    ΔU = windspeed(sc)
    
    # Ri_b = (g z (DSE_v,in - DSE_v,sfc)) / (DSE_v,sfc |u(z)|^2)
    numerator = 𝑔 * z * (DSEᵥ_in_val - DSEᵥ_sfc_val)
    denominator = DSEᵥ_sfc_val * ΔU^2
    
    if abs(denominator) < eps(FT)
        return FT(0)  # Neutral conditions
    end
    return numerator / denominator
end

"""
    bulk_richardson_number(param_set, sc, X★, scheme)

Compute the bulk Richardson number from similarity scales.

The bulk Richardson number is defined as:
    Ri_b = ζ F_c(ζ) F_m(ζ)^(-2)

where ζ = Δz / L★ is the Monin-Obukhov stability parameter,
F_c and F_m are the stability correction functions for scalar and momentum transport.
"""
function bulk_richardson_number(param_set, sc, X★, scheme)
    uf = SFP.uf_params(param_set)
    𝜅 = SFP.von_karman_const(param_set)
    FT = eltype(X★.L★)
    
    # Compute stability parameter
    ζ = Δz(sc) / X★.L★
    
    # Compute stability correction functions
    F_c = compute_Fₘₕ(sc, uf, ζ, X★.𝓁θ, UF.HeatTransport())
    F_m = compute_Fₘₕ(sc, uf, ζ, X★.𝓁u, UF.MomentumTransport())
    
    # Ri_b = ζ F_c(ζ) F_m(ζ)^(-2)
    # Avoid division by zero
    if abs(F_m) < eps(FT)
        return FT(Inf) * sign(ζ)
    end
    return ζ * F_c / (F_m^2)
end

"""
    bulk_richardson_number_from_zeta(param_set, sc, ζ, 𝓁u, 𝓁θ, scheme)

Compute the bulk Richardson number from stability parameter ζ and roughness lengths.

This is used for root finding where we solve for ζ.
"""
function bulk_richardson_number_from_zeta(param_set, sc, ζ, 𝓁u, 𝓁θ, scheme)
    uf = SFP.uf_params(param_set)
    FT = eltype(ζ)
    
    # Compute stability correction functions
    F_c = compute_Fₘₕ(sc, uf, ζ, 𝓁θ, UF.HeatTransport())
    F_m = compute_Fₘₕ(sc, uf, ζ, 𝓁u, UF.MomentumTransport())
    
    # Ri_b = ζ F_c(ζ) F_m(ζ)^(-2)
    # Avoid division by zero
    if abs(F_m) < eps(FT)
        return FT(Inf) * sign(ζ)
    end
    return ζ * F_c / (F_m^2)
end

"""
    obukhov_iteration(X★, sc, scheme, param_set, tol, maxiter, solver_method)

Iterative solver for Monin-Obukhov similarity solution using fixed point iteration
with convergence based on bulk Richardson number.

## Arguments
- `X★`: Initial guess for similarity scales (NamedTuple)
- `sc`: Surface conditions container
- `scheme`: Discretization scheme
- `param_set`: Parameter set
- `tol`: Convergence tolerance
- `maxiter`: Maximum number of iterations
- `solver_method`: Solver method (FixedPointIteration, BrentsMethod, or SecantMethod)

## Returns
Converged similarity scales (NamedTuple)
"""
function obukhov_iteration(X★,
    sc,
    scheme,
    param_set,
    tol,
    maxiter = 20,
    solver_method::SolverMethod = FixedPointIteration(),
)
    return obukhov_iteration_fixed_point(X★, sc, scheme, param_set, tol, maxiter)
end

function obukhov_iteration_fixed_point(X★,
    sc,
    scheme,
    param_set,
    tol,
    maxiter,
)
    FT = eltype(X★)
    qₛ = surface_specific_humidity(param_set, sc)
    
    # Define fixed point function
    function fixed_point_func(X★_in)
        return iterate_interface_fluxes(sc,
            qₛ,
            X★_in,
            ts_in(sc),
            ts_sfc(sc),
            scheme,
            param_set)
    end
    
    # Define convergence check using bulk Richardson number
    function convergence_check(X★_prev, X★_curr)
        Ri_b_prev = bulk_richardson_number(param_set, sc, X★_prev, scheme)
        Ri_b_curr = bulk_richardson_number(param_set, sc, X★_curr, scheme)
        return abs(Ri_b_curr - Ri_b_prev) ≤ tol
    end
    
    # Use fixed point iteration with custom convergence check
    X★_prev = X★
    for iter in 1:maxiter
        X★_curr = fixed_point_func(X★_prev)
        if convergence_check(X★_prev, X★_curr)
            return X★_curr
        end
        X★_prev = X★_curr
    end
    return X★_prev
end

"""
    reconstruct_X★_from_zeta(param_set, sc, ζ, scheme, X★_guess)

Reconstruct similarity scales X★ from stability parameter ζ.
Uses an iterative approach to handle roughness length dependence on u★.
"""
function reconstruct_X★_from_zeta(param_set, sc, ζ, scheme, X★_guess)
    uf = SFP.uf_params(param_set)
    𝜅 = SFP.von_karman_const(param_set)
    FT = eltype(ζ)
    
    # Compute L★ from ζ
    L★ = Δz(sc) / ζ
    
    # Initial guess for u★ (needed for roughness lengths)
    u★ = X★_guess.u★
    ΔU = windspeed(sc)
    
    # Iterate to get consistent u★ and roughness lengths
    for _ in 1:5  # Small number of iterations for consistency
        𝓁u = compute_z0(u★, param_set, sc, sc.roughness_model, UF.MomentumTransport())
        𝓁θ = compute_z0(u★, param_set, sc, sc.roughness_model, UF.HeatTransport())
        𝓁q = compute_z0(u★, param_set, sc, sc.roughness_model, UF.HeatTransport())
        
        # Compute χ coefficients
        F_m = compute_Fₘₕ(sc, uf, ζ, 𝓁u, UF.MomentumTransport())
        χu = 𝜅 / F_m
        u★ = χu * ΔU
    end
    
    # Final computation of all scales
    𝓁u = compute_z0(u★, param_set, sc, sc.roughness_model, UF.MomentumTransport())
    𝓁θ = compute_z0(u★, param_set, sc, sc.roughness_model, UF.HeatTransport())
    𝓁q = compute_z0(u★, param_set, sc, sc.roughness_model, UF.HeatTransport())
    
    F_m = compute_Fₘₕ(sc, uf, ζ, 𝓁u, UF.MomentumTransport())
    F_c = compute_Fₘₕ(sc, uf, ζ, 𝓁θ, UF.HeatTransport())
    F_q = compute_Fₘₕ(sc, uf, ζ, 𝓁q, UF.HeatTransport())
    
    χu = 𝜅 / F_m
    χDSEᵥ = 𝜅 / F_c
    χq = 𝜅 / F_q
    χθᵥ = 𝜅 / F_c
    
    u★ = χu * ΔU
    DSEᵥ★ = χDSEᵥ * ΔDSEᵥ(param_set, sc)
    q★ = χq * Δqt(param_set, sc)
    θᵥ★ = χθᵥ * Δθᵥ(param_set, sc)
    
    return (; u★, DSEᵥ★, q★, L★, θᵥ★, 𝓁u, 𝓁θ, 𝓁q)
end

# Method dispatch for different solver methods
# Support RootSolvers.jl methods directly
function obukhov_iteration(X★,
    sc,
    scheme,
    param_set,
    tol,
    maxiter,
    solver_method::RS.BrentsMethod,
)
    return obukhov_iteration_rootsolver(X★, sc, scheme, param_set, tol, maxiter, solver_method)
end

function obukhov_iteration(X★,
    sc,
    scheme,
    param_set,
    tol,
    maxiter,
    solver_method::RS.SecantMethod,
)
    return obukhov_iteration_rootsolver(X★, sc, scheme, param_set, tol, maxiter, solver_method)
end

# Backwards compatibility with marker types
function obukhov_iteration(X★,
    sc,
    scheme,
    param_set,
    tol,
    maxiter,
    solver_method::BrentsMethod,
)
    # Create RootSolvers method with auto-computed brackets
    return obukhov_iteration_rootsolver(X★, sc, scheme, param_set, tol, maxiter, nothing, :brent)
end

function obukhov_iteration(X★,
    sc,
    scheme,
    param_set,
    tol,
    maxiter,
    solver_method::SecantMethod,
)
    # Create RootSolvers method with auto-computed initial guesses
    return obukhov_iteration_rootsolver(X★, sc, scheme, param_set, tol, maxiter, nothing, :secant)
end

"""
    obukhov_iteration_rootsolver(X★, sc, scheme, param_set, tol, maxiter, rootsolver_method)

Use RootSolvers.jl to solve for stability parameter ζ by finding the root of the
Richardson number residual function.

The residual function is: f(ζ) = Ri_b(ζ) - Ri_b_target, where Ri_b(ζ) is computed
from the stability correction functions and Ri_b_target is computed from the input
state variables.

## Arguments
- `rootsolver_method`: Either a `RootSolvers.BrentsMethod` or `RootSolvers.SecantMethod` instance,
  or `nothing` with `method_type` symbol for auto-computed brackets.
"""
function obukhov_iteration_rootsolver(X★,
    sc,
    scheme,
    param_set,
    tol,
    maxiter,
    rootsolver_method::Union{RS.BrentsMethod, RS.SecantMethod},
)
    return obukhov_iteration_rootsolver(X★, sc, scheme, param_set, tol, maxiter, rootsolver_method, nothing)
end

function obukhov_iteration_rootsolver(X★,
    sc,
    scheme,
    param_set,
    tol,
    maxiter,
    rootsolver_method::Union{RS.BrentsMethod, RS.SecantMethod, Nothing},
    method_type::Union{Symbol, Nothing},
)
    FT = eltype(X★)
    
    # Compute target Richardson number from input state
    Ri_b_target = bulk_richardson_number_target(param_set, sc)
    
    # Initial guess for ζ
    ζ₀ = Δz(sc) / X★.L★
    
    # Get initial roughness lengths for first iteration
    u★_init = X★.u★
    𝓁u_init = compute_z0(u★_init, param_set, sc, sc.roughness_model, UF.MomentumTransport())
    𝓁θ_init = compute_z0(u★_init, param_set, sc, sc.roughness_model, UF.HeatTransport())
    
    uf = SFP.uf_params(param_set)
    𝜅 = SFP.von_karman_const(param_set)
    ΔU = windspeed(sc)
    
    # Define residual function: f(ζ) = Ri_b(ζ) - Ri_b_target
    # For each ζ, we need to compute Ri_b using the stability correction functions
    # Note: roughness lengths depend on u★ which depends on ζ, so we iterate
    # to get consistent values within the residual evaluation
    function residual(ζ)
        # Iterate to get consistent u★ and roughness lengths for this ζ
        u★ = u★_init
        for _ in 1:3  # Small number of iterations for consistency
            𝓁u = compute_z0(u★, param_set, sc, sc.roughness_model, UF.MomentumTransport())
            𝓁θ = compute_z0(u★, param_set, sc, sc.roughness_model, UF.HeatTransport())
            
            # Compute F_m and update u★
            F_m = compute_Fₘₕ(sc, uf, ζ, 𝓁u, UF.MomentumTransport())
            if abs(F_m) < eps(FT)
                return FT(Inf) * sign(ζ)  # Avoid division by zero
            end
            χu = 𝜅 / F_m
            u★ = χu * ΔU
        end
        
        # Final computation with consistent roughness lengths
        𝓁u = compute_z0(u★, param_set, sc, sc.roughness_model, UF.MomentumTransport())
        𝓁θ = compute_z0(u★, param_set, sc, sc.roughness_model, UF.HeatTransport())
        
        # Compute Ri_b from ζ
        Ri_b_ζ = bulk_richardson_number_from_zeta(param_set, sc, ζ, 𝓁u, 𝓁θ, scheme)
        return Ri_b_ζ - Ri_b_target
    end
    
    # Use RootSolvers to find ζ
    try
        # Use provided method or create one based on method_type
        if rootsolver_method !== nothing
            method = rootsolver_method
        else
            # Set up initial bracket/guesses based on stability
            # For stable conditions: ζ > 0, for unstable: ζ < 0
            if ζ₀ >= 0
                # Stable conditions: ζ > 0
                bracket_low = max(FT(1e-6), ζ₀ * FT(0.1))
                bracket_high = FT(10.0)
            else
                # Unstable conditions: ζ < 0
                bracket_low = FT(-10.0)
                bracket_high = min(FT(-1e-6), ζ₀ * FT(0.1))
            end
            
            if method_type === :brent
                # Brent's method requires a bracketing interval
                method = RS.BrentsMethod(bracket_low, bracket_high)
            elseif method_type === :secant
                # Secant method requires two initial guesses
                x0 = bracket_low
                x1 = bracket_high
                method = RS.SecantMethod(x0, x1)
            else
                error("Unknown method type: $method_type")
            end
        end
        
        # Find the root
        sol = RS.find_zero(residual, method; x_tol = tol, max_iters = maxiter)
        
        # Check convergence
        if RS.converged(sol)
            ζ_solved = sol.x
            # Reconstruct X★ from solved ζ
            return reconstruct_X★_from_zeta(param_set, sc, ζ_solved, scheme, X★)
        else
            # Fall back to fixed point iteration if root finding fails
            return obukhov_iteration_fixed_point(X★, sc, scheme, param_set, tol, maxiter)
        end
    catch e
        # Fall back to fixed point iteration if root finding fails
        return obukhov_iteration_fixed_point(X★, sc, scheme, param_set, tol, maxiter)
    end
end

# For backwards compatibility with package extensions
if !isdefined(Base, :get_extension)
    include(joinpath("..", "ext", "CreateParametersExt.jl"))
end

end # SurfaceFluxes module
