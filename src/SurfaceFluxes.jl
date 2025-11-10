"""
    SurfaceFluxes

Surface flux parameterization module implementing Monin-Obukhov similarity theory
for computing surface exchange coefficients and fluxes in atmospheric boundary layer models.

## Features

  - Monin-Obukhov similarity theory implementation
  - Support for multiple discretization schemes (finite difference and finite volume)
  - Multiple surface condition types (prescribed fluxes, coefficients, or values only)
  - Roughness parameterization options (scalar and Charnock)
  - Universal function support for stable/unstable boundary layers

## Main Functions

  - [`surface_conditions`](@ref): Main interface for computing surface conditions
  - Computes Monin-Obukhov length
  - Computes sensible and latent heat fluxes
  - Computes momentum fluxes
  - Computes exchange coefficients (Cd, Ch)
  - Computes friction velocity (u★)

## Discretization Schemes

  - `PointValueScheme`: Finite difference scheme (Byun 1990)
  - `LayerAverageScheme`: Finite volume scheme (Nishizawa 2018)

## References
 - [Nishizawa2018](@cite): Finite volume discretization scheme
 - [Byun1990](@cite): Finite difference discretization scheme

"""
module SurfaceFluxes


import RootSolvers
const RS = RootSolvers

using DocStringExtensions
const DSE = DocStringExtensions

import Thermodynamics
const TD = Thermodynamics

include("UniversalFunctions.jl")
import .UniversalFunctions
const UF = UniversalFunctions

include("Parameters.jl")
import .Parameters

const SFP = Parameters
const APS = SFP.AbstractSurfaceFluxesParameters

"""Abstract type for discretization schemes."""
abstract type SolverScheme end

"""Finite volume scheme (Nishizawa 2018)."""
struct LayerAverageScheme <: SolverScheme end

"""Finite difference scheme (Byun 1990)."""
struct PointValueScheme <: SolverScheme end

"""
    SurfaceFluxConditions

Surface flux conditions, returned from `surface_conditions`.

# Fields

$(DSE.FIELDS)
"""
struct SurfaceFluxConditions{FT <: Real}
    L_MO::FT
    shf::FT
    lhf::FT
    buoy_flux::FT
    ρτxz::FT
    ρτyz::FT
    ustar::FT
    Cd::FT
    Ch::FT
    evaporation::FT
end

SurfaceFluxConditions(L_MO, shf, lhf, buoy_flux, ρτxz, ρτyz, ustar, Cd, Ch, E) =
    SurfaceFluxConditions(
        promote(L_MO, shf, lhf, buoy_flux, ρτxz, ρτyz, ustar, Cd, Ch, E)...,
    )

function Base.show(io::IO, sfc::SurfaceFluxConditions)
    println(io, "----------------------- SurfaceFluxConditions")
    println(io, "L_MO                   = ", sfc.L_MO)
    println(io, "Sensible Heat Flux     = ", sfc.shf)
    println(io, "Latent Heat Flux       = ", sfc.lhf)
    println(io, "Buoyancy Flux          = ", sfc.buoy_flux)
    println(io, "Friction velocity u⋆   = ", sfc.ustar)
    println(io, "C_drag                 = ", sfc.Cd)
    println(io, "C_heat                 = ", sfc.Ch)
    println(io, "evaporation            = ", sfc.evaporation)
    println(io, "-----------------------")
end

"""Abstract type for roughness models."""
abstract type RoughnessModel end

"""Constant momentum and scalar roughness lengths."""
struct ScalarRoughness <: RoughnessModel end

"""Charnock's relation for momentum roughness."""
struct CharnockRoughness <: RoughnessModel end

"""
    MOSTSolutionState

Internal mutable struct used to capture intermediate values during Brent's method iteration.
"""
mutable struct MOSTSolutionState{FT}
    ustar::FT
    z0m::FT
end

"""
   StateValues

Input container for state variables at either first / interior nodes.

# Fields

$(DSE.FIELDS)
"""
struct StateValues{FT <: Real, A, TS <: TD.ThermodynamicState}
    z::FT
    u::A
    ts::TS
end

"""Abstract type for surface condition containers."""
abstract type AbstractSurfaceConditions{
    FT <: Real,
    SVA <: StateValues,
    SVB <: StateValues,
    RM <: RoughnessModel,
} end


"""
    Fluxes

Input container for state variables, latent and sensible heat fluxes roughness lengths,
initial obukhov length and gustiness.

# Fields

$(DSE.FIELDS)
"""
struct Fluxes{FT, SVA, SVB, RM} <: AbstractSurfaceConditions{FT, SVA, SVB, RM}
    state_in::SVA
    state_sfc::SVB
    shf::FT
    lhf::FT
    z0m::FT
    z0b::FT
    gustiness::FT
    roughness_model::RM
end

function Fluxes(
    state_in::SVA,
    state_sfc::SVB,
    shf::FT,
    lhf::FT,
    z0m::FT,
    z0b::FT;
    gustiness::FT = FT(1),
    roughness_model::RM = ScalarRoughness(),
) where {SVA, SVB, FT, RM}
    return Fluxes{FT, SVA, SVB, RM}(
        state_in,
        state_sfc,
        shf,
        lhf,
        z0m,
        z0b,
        gustiness,
        roughness_model,
    )
end


"""
    FluxesAndFrictionVelocity

Input container, given surface state variables, latent and sensible heat fluxes,
and the friction velocity, roughness lengths,
initial obukhov length and gustiness.

# Fields

$(DSE.FIELDS)
"""
struct FluxesAndFrictionVelocity{FT, SVA, SVB, RM} <:
       AbstractSurfaceConditions{FT, SVA, SVB, RM}
    state_in::SVA
    state_sfc::SVB
    shf::FT
    lhf::FT
    ustar::FT
    z0m::FT
    z0b::FT
    gustiness::FT
    roughness_model::RM
end

function FluxesAndFrictionVelocity(
    state_in::SVA,
    state_sfc::SVB,
    shf::FT,
    lhf::FT,
    ustar::FT,
    z0m::FT,
    z0b::FT;
    gustiness::FT = FT(1),
    roughness_model::RM = ScalarRoughness(),
) where {SVA, SVB, FT, RM}
    return FluxesAndFrictionVelocity{FT, SVA, SVB, RM}(
        state_in,
        state_sfc,
        shf,
        lhf,
        ustar,
        z0m,
        z0b,
        gustiness,
        roughness_model,
    )
end

"""
    Coefficients

Input container, given surface state variables, and exchange coefficients,roughness lengths,
initial obukhov length and gustiness.

# Fields

$(DSE.FIELDS)
"""
struct Coefficients{FT, SVA, SVB, RM} <: AbstractSurfaceConditions{FT, SVA, SVB, RM}
    state_in::SVA
    state_sfc::SVB
    Cd::FT
    Ch::FT
    gustiness::FT
    beta::FT
    roughness_model::RM
end

function Coefficients(
    state_in::SVA,
    state_sfc::SVB,
    Cd::FT,
    Ch::FT;
    gustiness::FT = FT(1),
    beta::FT = FT(1),
    roughness_model::RM = ScalarRoughness(),
) where {FT, SVA, SVB, RM}
    return Coefficients{FT, SVA, SVB, RM}(
        state_in,
        state_sfc,
        Cd,
        Ch,
        gustiness,
        beta,
        roughness_model,
    )
end


"""
    ValuesOnly

Input container, given only surface state variables, roughness lengths,
initial obukhov length and gustiness.

# Fields

$(DSE.FIELDS)
"""
struct ValuesOnly{FT, SVA, SVB, RM} <: AbstractSurfaceConditions{FT, SVA, SVB, RM}
    state_in::SVA
    state_sfc::SVB
    z0m::FT
    z0b::FT
    gustiness::FT
    beta::FT
    roughness_model::RM
end

function ValuesOnly(
    state_in::SVA,
    state_sfc::SVB,
    z0m::FT,
    z0b::FT;
    gustiness::FT = FT(1),
    beta::FT = FT(1),
    roughness_model::RM = ScalarRoughness(),
) where {SVA, SVB, FT, RM}
    return ValuesOnly{FT, SVA, SVB, RM}(
        state_in,
        state_sfc,
        z0m,
        z0b,
        gustiness,
        beta,
        roughness_model,
    )
end

"""Thermodynamic state at interior level."""
ts_in(sc::AbstractSurfaceConditions) = sc.state_in.ts

"""Thermodynamic state at surface."""
ts_sfc(sc::AbstractSurfaceConditions) = sc.state_sfc.ts

"""Height at interior level (m)."""
z_in(sc::AbstractSurfaceConditions) = sc.state_in.z

"""Height at surface level (m)."""
z_sfc(sc::AbstractSurfaceConditions) = sc.state_sfc.z

"""Vertical displacement between interior and surface (m)."""
Δz(sc::AbstractSurfaceConditions) = z_in(sc) - z_sfc(sc)

z0(sc, transport) = z0(sc, transport, sc.roughness_model)

z0(sc::AbstractSurfaceConditions, ::UF.MomentumTransport, ::ScalarRoughness) = sc.z0m
z0(sc::AbstractSurfaceConditions, ::UF.HeatTransport, ::ScalarRoughness) = sc.z0b
z0(sc::AbstractSurfaceConditions, ::UF.HeatTransport, ::CharnockRoughness) = sc.z0b

z0(sc::Coefficients, ::UF.MomentumTransport, ::ScalarRoughness) = nothing
z0(sc::Coefficients, ::UF.HeatTransport, ::ScalarRoughness) = nothing

Δu1(sc::AbstractSurfaceConditions) = sc.state_in.u[1] - sc.state_sfc.u[1]
Δu2(sc::AbstractSurfaceConditions) = sc.state_in.u[2] - sc.state_sfc.u[2]

"""Total specific humidity at interior level."""
qt_in(param_set::APS, sc::AbstractSurfaceConditions) =
    TD.total_specific_humidity(SFP.thermodynamics_params(param_set), ts_in(sc))

"""Total specific humidity at surface level."""
qt_sfc(param_set::APS, sc::AbstractSurfaceConditions) =
    TD.total_specific_humidity(SFP.thermodynamics_params(param_set), ts_sfc(sc))

"""Specific humidity difference between interior and surface levels."""
Δqt(param_set::APS, sc::AbstractSurfaceConditions) =
    qt_in(param_set, sc) - qt_sfc(param_set, sc)

"""Virtual dry static energy at interior level."""
DSEᵥ_in(param_set::APS, sc::AbstractSurfaceConditions) =
    TD.virtual_dry_static_energy(SFP.thermodynamics_params(param_set),
        ts_in(sc),
        SFP.grav(param_set) * z_in(sc))

"""Virtual dry static energy at surface level."""
DSEᵥ_sfc(param_set::APS, sc::AbstractSurfaceConditions) =
    TD.virtual_dry_static_energy(
        SFP.thermodynamics_params(param_set),
        ts_sfc(sc),
        SFP.grav(param_set) * z_sfc(sc),
    )

"""Virtual dry static energy difference between interior and surface levels."""
ΔDSEᵥ(param_set, sc) = DSEᵥ_in(param_set, sc) - DSEᵥ_sfc(param_set, sc)

"""Wind velocity vector at interior level."""
u_in(sc::AbstractSurfaceConditions) = sc.state_in.u

"""Wind velocity vector at surface level."""
u_sfc(sc::AbstractSurfaceConditions) = sc.state_sfc.u

"""Wind speed magnitude with gustiness."""
function windspeed(sc::AbstractSurfaceConditions)
    return max(hypot(Δu1(sc), Δu2(sc)), sc.gustiness)
end

"""
    surface_conditions(
        param_set::AbstractSurfaceFluxesParameters,
        sc::SurfaceFluxes.AbstractSurfaceConditions,
        scheme::SurfaceFluxes.SolverScheme = PointValueScheme();
        tol_neutral = SFP.cp_d(param_set) / 100,
        tol::RS.AbstractTolerance = RS.RelativeSolutionTolerance(FT(0.01)),
        maxiter::Int = 10,
        soltype::RS.SolutionType = RS.CompactSolution(),
    )

Main interface for computing surface conditions using Monin-Obukhov similarity theory.

# Arguments
- `param_set`: Parameter set containing physical and thermodynamic constants
- `sc`: Surface conditions container (Fluxes, FluxesAndFrictionVelocity, Coefficients, or ValuesOnly)
- `scheme`: Discretization scheme (PointValueScheme or LayerAverageScheme), default: PointValueScheme()
- `tol_neutral`: Tolerance for neutral boundary layer detection, default: cp_d/100
- `tol`: Root solver tolerance, default: RelativeSolutionTolerance(0.01)
- `maxiter`: Maximum iterations for root solver, default: 10
- `soltype`: Solution type from RootSolvers, default: CompactSolution()

# Returns
Returns `SurfaceFluxConditions` struct containing:
  - `L_MO`: Monin-Obukhov lengthscale (m)
  - `shf`: Sensible heat flux (W/m²)
  - `lhf`: Latent heat flux (W/m²)
  - `ρτxz`: Momentum flux eastward component (kg/m/s²)
  - `ρτyz`: Momentum flux northward component (kg/m/s²)
  - `ustar`: Friction velocity (m/s)
  - `Cd`: Momentum exchange coefficient
  - `Ch`: Thermal exchange coefficient
  - `evaporation`: Evaporation rate (kg/m²/s)
"""
function surface_conditions(
    param_set::APS{FT},
    sc::AbstractSurfaceConditions,
    scheme::SolverScheme = PointValueScheme();
    tol_neutral = SFP.cp_d(param_set) / 100,
    tol::RS.AbstractTolerance = RS.RelativeSolutionTolerance(FT(0.01)),
    maxiter::Int = 10,
    soltype::RS.SolutionType = RS.CompactSolution(),
) where {FT}
    L_MO, z0m, ustar = obukhov_similarity_solution(
        param_set,
        sc,
        scheme,
        tol,
        tol_neutral,
        maxiter,
        soltype,
    )
    z0b = z0(sc, UF.HeatTransport())

    Cd = momentum_exchange_coefficient(
        param_set,
        L_MO,
        sc,
        scheme,
        tol_neutral,
        z0m,
        z0b,
    )
    Ch = heat_exchange_coefficient(
        param_set,
        L_MO,
        sc,
        scheme,
        tol_neutral,
        z0m,
        z0b,
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
    sfc_results = SurfaceFluxConditions(
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

"""
    obukhov_similarity_solution(sfc::SurfaceFluxConditions)

    obukhov_similarity_solution(
        param_set::AbstractSurfaceFluxesParameters,
        sc::AbstractSurfaceConditions,
        scheme,
        tol,
        tol_neutral,
        maxiter,
        soltype,
    )

Iteratively solves for the Monin-Obukhov lengthscale, aerodynamic roughness
and friction velocity u★.

# Arguments
- `param_set`: Parameter set with physical constants
- `sc`: Surface conditions container
- `uft`: Universal function type
- `scheme`: Discretization scheme
- `tol`: Root solver tolerance
- `tol_neutral`: Tolerance for neutral layer detection
- `maxiter`: Maximum iterations for root solver
- `soltype`: Root solver solution type

# Returns
Returns `(L_MO, z0m, ustar)` tuple where:
  - `L_MO`: Monin-Obukhov length (m)
  - `z0m`: Momentum roughness length (m)
  - `ustar`: Friction velocity (m/s)

# Methods
For `FluxesAndFrictionVelocity`: Direct computation, no iteration needed.
For `Coefficients`: Computed from exchange coefficients, no iteration needed.
For `Fluxes` or `ValuesOnly`: Iterative solution using root solver.
"""
function obukhov_similarity_solution end

function obukhov_similarity_solution(sfc::SurfaceFluxConditions)
    return sfc.L_MO
end

"""Ensure non-zero value (prevents division by zero)."""
function non_zero(v::FT) where {FT}
    sign_of_v = v == 0 ? 1 : sign(v)
    return abs(v) < eps(FT) ? eps(FT) * sign_of_v : v
end

"""Compute bulk Richardson number given the atmospheric state properties at the 
interior and surface levels."""
function compute_richardson_number(
    sc::AbstractSurfaceConditions,
    DSEᵥ_in,
    DSEᵥ_sfc,
    grav,
)
    return (grav * Δz(sc) * (DSEᵥ_in - DSEᵥ_sfc)) /
           (DSEᵥ_in * (windspeed(sc))^2)
end

"""
    compute_charnock_roughness(param_set, u★)

Compute momentum roughness length using Charnock's relation (Charnock 1955).

# Arguments
- `param_set`: Parameter set with physical constants
- `u★`: Friction velocity (m/s)

# Returns
Momentum roughness length (m).

Computed as: z0m = α_charnock * u★² / g, where α_charnock = 0.011.
"""
function compute_charnock_roughness(param_set, u★)
    FT = eltype(u★)
    α_charnock = FT(0.011)
    return α_charnock * u★^2 / SFP.grav(param_set)
end
function compute_Ri_b(param_set, sc::AbstractSurfaceConditions,  scheme, ζ, ::ScalarRoughness, solution_state)
    thermo_params = SFP.thermodynamics_params(param_set)
    ufparams = SFP.uf_params(param_set)
    𝓏0m = z0(sc, UF.MomentumTransport())
    𝓏0b = z0(sc, UF.HeatTransport())
    ufₛ = SFP.uf_params(param_set)
    u★ = compute_ustar(param_set, Δz(sc) / ζ, sc,  scheme, 𝓏0m, 𝓏0b)
    F_m = compute_Fₘₕ(sc, ufₛ, ζ, 𝓏0m, UF.MomentumTransport())
    F_h = compute_Fₘₕ(sc, ufₛ, ζ, 𝓏0b, UF.HeatTransport())
    return (ζ * F_h / F_m^2, u★, 𝓏0m)
end
function compute_Ri_b(param_set, sc::AbstractSurfaceConditions,  scheme, ζ, ::CharnockRoughness, solution_state)
    thermo_params = SFP.thermodynamics_params(param_set)
    ufparams = SFP.uf_params(param_set)
    u★ = solution_state.ustar
    𝓏0m = compute_charnock_roughness(param_set, u★)
    ufₛ = SFP.uf_params(param_set)
    F_m = compute_Fₘₕ(sc, ufₛ, ζ, 𝓏0m, UF.MomentumTransport())
    F_h = compute_Fₘₕ(sc, ufₛ, ζ, sc.z0b, UF.HeatTransport())
    return (ζ * F_h / F_m^2, u★, 𝓏0m)
end

"""
    compute_Fₘₕ(sc::AbstractSurfaceConditions, ufₛ, ζ, 𝓏0, transport)

Compute integrated universal function F_m or F_h.
    log(z(sc)/z₀) - ψ(z/L) + ψ(z₀/L)

# Arguments
- `sc`: Surface conditions container, with model level given by z(sc)
- `ufₛ`: Universal function
- `ζ`: Stability parameter
- `𝓏0`: Aerodynamic roughness length
-
`transport`: Transport type
"""
function compute_Fₘₕ(sc::AbstractSurfaceConditions, ufₛ, ζ, 𝓏0, transport)
    ψ = UF.psi(ufₛ, ζ, transport)
    ψ₀ = UF.psi(ufₛ, 𝓏0 * ζ / Δz(sc), transport)
    return log(Δz(sc) / 𝓏0) - ψ + ψ₀
end

function obukhov_similarity_solution(
    param_set::APS{FT},
    sc::Union{Fluxes, ValuesOnly},
    scheme,
    tol,
    tol_neutral,
    maxiter,
    soltype,
) where {FT}
    thermo_params = SFP.thermodynamics_params(param_set)
    ufparams = SFP.uf_params(param_set)
    grav = SFP.grav(param_set)
    𝜅 = SFP.von_karman_const(param_set)

    # Capture the final values during iteration
    solution_state = MOSTSolutionState(FT(0.1), FT(0.001))
    function root_ζ(ζ)
        f1 = compute_richardson_number(sc,
            DSEᵥ_in(param_set, sc),
            DSEᵥ_sfc(param_set, sc),
            grav)
        f2, u★, 𝓏0 = compute_Ri_b(param_set, sc,  scheme, ζ,
            sc.roughness_model,
            solution_state)
        # Capture the values from this iteration
        solution_state.ustar = non_zero(u★)
        solution_state.z0m = non_zero(𝓏0)
        return f1 - f2
    end
    function root_u★(u★)
        f1 = windspeed(sc)
        f2 = non_zero(u★) / 𝜅 * 
             log(Δz(sc) / non_zero(compute_charnock_roughness(param_set, u★)))
        return f1 - f2
    end

    if ΔDSEᵥ(param_set, sc) >= FT(0)
        # Iterative Solution where ζ > 0 (Stable BL)
        sol = RS.find_zero(
            root_ζ,
            RS.BrentsMethod(FT(0), FT(1e6)),
            soltype,
            tol,
            maxiter,
        )
        ζ = sol.root
        L_MO = Δz(sc) / non_zero(ζ)
        return (non_zero(L_MO), solution_state.z0m, solution_state.ustar)
    else
        # Iterative Solution where ζ < 0 (Unstable BL)
        sol = RS.find_zero(
            root_ζ,
            RS.BrentsMethod(FT(-1e6), FT(0)),
            soltype,
            tol,
            maxiter,
        )
        ζ = sol.root
        L_MO = Δz(sc) / non_zero(ζ)
        return (non_zero(L_MO), solution_state.z0m, solution_state.ustar)
    end
end

function obukhov_similarity_solution(
    param_set,
    sc::FluxesAndFrictionVelocity,
    scheme,
    args...,
)
    L_MO = -sc.ustar^3 / SFP.von_karman_const(param_set) /
           non_zero(compute_buoyancy_flux(param_set, sc, scheme))
    ustar_val = sc.ustar
    z0m_val = z0(sc, UF.MomentumTransport())
    return (L_MO, z0m_val, ustar_val)
end

function obukhov_similarity_solution(
    param_set,
    sc::Coefficients,
    scheme,
    args...,
)
    lhf = latent_heat_flux(param_set, sc.Ch, sc, scheme)
    shf = sensible_heat_flux(param_set, sc.Ch, sc, scheme)
    ustar_val = sqrt(sc.Cd) * windspeed(sc)
    buoyancy_flux = compute_buoyancy_flux(
        param_set,
        shf,
        lhf,
        ts_in(sc),
        ts_sfc(sc),
        scheme,
    )
    L_MO = -ustar_val^3 / SFP.von_karman_const(param_set) / non_zero(buoyancy_flux)
    z0m_val = z0(sc, UF.MomentumTransport())
    return (L_MO, z0m_val, ustar_val)
end

"""
    compute_buoyancy_flux(param_set, shf, lhf, ts_in, ts_sfc, scheme)

Compute buoyancy flux from sensible and latent heat fluxes.

# Arguments
- `param_set`: Parameter set with physical constants
- `shf`: Sensible heat flux (W/m²)
- `lhf`: Latent heat flux (W/m²)
- `ts_in`: Thermodynamic state at interior level
- `ts_sfc`: Thermodynamic state at surface
- `scheme`: Discretization scheme

# Returns
Buoyancy flux (m²/s³).

Computed as: B = (g/ρ_sfc) * (shf/(cp_m * T_in) + (ε_vd - 1) * lhf/L_v)
"""
function compute_buoyancy_flux(param_set, shf, lhf, ts_in, ts_sfc, scheme)
    thermo_params = SFP.thermodynamics_params(param_set)
    grav = SFP.grav(param_set)
    ε_vd = SFP.Rv_over_Rd(param_set)
    cp_m = TD.cp_m(thermo_params, ts_in)
    L_v = TD.latent_heat_vapor(thermo_params, ts_in)
    ρ_sfc = TD.air_density(thermo_params, ts_sfc)
    T_in = TD.air_temperature(thermo_params, ts_in)
    return grav / ρ_sfc * (shf / cp_m / T_in + (ε_vd - 1) * lhf / L_v)
end

function compute_buoyancy_flux(
    param_set,
    sc::Union{FluxesAndFrictionVelocity, Fluxes},
    scheme,
)
    return compute_buoyancy_flux(
        param_set,
        sc.shf,
        sc.lhf,
        ts_in(sc),
        ts_sfc(sc),
        scheme,
    )
end

"""
    compute_ustar(
        param_set::AbstractSurfaceFluxesParameters,
        L_MO,
        sc::AbstractSurfaceConditions,
        
        scheme,
        args...
    )

Compute friction velocity (u★).

# Arguments
- `param_set`: Parameter set with physical constants
- `L_MO`: Monin-Obukhov length (m)
- `sc`: Surface conditions container
- `uft`: Universal function type
- `scheme`: Discretization scheme
- `args...`: Optional arguments (e.g., z0m, z0b)

# Returns
Friction velocity (m/s).

# Methods
- `FluxesAndFrictionVelocity`: Returns known friction velocity from `sc`.
- `Fluxes` or `ValuesOnly`: Computed from Monin-Obukhov length and universal functions.
- `Coefficients`: Computed as `√(Cd) * windspeed(sc)`.
"""
function compute_ustar end

compute_ustar(param_set, L_MO, sc::FluxesAndFrictionVelocity,  scheme, args...) =
    sc.ustar

function compute_ustar(param_set, L_MO, sc::Fluxes,  scheme, z0m = nothing, z0b = nothing)
    z0m_val = z0m === nothing ? z0(sc, UF.MomentumTransport()) : z0m
    z0b_val = z0b === nothing ? z0(sc, UF.HeatTransport()) : z0b
    windspeed(sc) * compute_physical_scale_coeff(
        param_set,
        sc,
        L_MO,
        UF.MomentumTransport(),
        scheme,
        z0m_val,
        z0b_val,
    )
end

compute_ustar(param_set, L_MO, sc::Coefficients,  scheme, args...) =
    sqrt(sc.Cd) * (windspeed(sc))

function compute_ustar(param_set, L_MO, sc::ValuesOnly,  scheme, z0m = nothing, z0b = nothing)
    z0m_val = z0m === nothing ? z0(sc, UF.MomentumTransport()) : z0m
    z0b_val = z0b === nothing ? z0(sc, UF.HeatTransport()) : z0b
    windspeed(sc) * compute_physical_scale_coeff(
        param_set,
        sc,
        L_MO,
        UF.MomentumTransport(),
        scheme,
        z0m_val,
        z0b_val,
    )
end

"""
    momentum_exchange_coefficient(param_set, L_MO, sc,  scheme, tol_neutral, z0m, z0b)

Compute momentum exchange coefficient (Cd).

# Arguments
- `param_set`: Parameter set with physical constants
- `L_MO`: Monin-Obukhov length (m)
- `sc`: Surface conditions container
- `uft`: Universal function type
- `scheme`: Discretization scheme
- `tol_neutral`: Tolerance for neutral layer detection
- `z0m`: Momentum roughness length (m)
- `z0b`: Scalar roughness length (m)

# Returns
Momentum exchange coefficient Cd (dimensionless).

For neutral conditions: Cd = (κ / log(Δz/z0m))².
For stable/unstable: Cd = u★² / windspeed².
"""
function momentum_exchange_coefficient(
    param_set,
    L_MO,
    sc::Union{Fluxes, ValuesOnly, FluxesAndFrictionVelocity},
    scheme,
    tol_neutral,
    z0m,
    z0b,
)
    thermo_params = SFP.thermodynamics_params(param_set)
    𝜅 = SFP.von_karman_const(param_set)
    grav = SFP.grav(param_set)
    if abs(ΔDSEᵥ(param_set, sc)) <= tol_neutral
        Cd = (𝜅 / log(Δz(sc) / z0m))^2
    else
        ustar = compute_ustar(param_set, L_MO, sc,  scheme, z0m, z0b)
        Cd = ustar^2 / windspeed(sc)^2
    end
    return Cd
end

"""
    momentum_exchange_coefficient(param_set, L_MO, sc,  scheme, tol_neutral, z0m, z0b)

Return Cd from Coefficients container.
"""
function momentum_exchange_coefficient(
    param_set,
    L_MO,
    sc::Coefficients,
    scheme,
    tol_neutral,
    z0m,
    z0b,
)
    return sc.Cd
end

"""
    heat_exchange_coefficient(param_set, L_MO, sc,  scheme, tol_neutral, z0m, z0b)

Compute heat exchange coefficient (Ch).

# Arguments
- `param_set`: Parameter set with physical constants
- `L_MO`: Monin-Obukhov length (m)
- `sc`: Surface conditions container
- `uft`: Universal function type
- `scheme`: Discretization scheme
- `tol_neutral`: Tolerance for neutral layer detection
- `z0m`: Momentum roughness length (m)
- `z0b`: Scalar roughness length (m)

# Returns
Heat exchange coefficient Ch (dimensionless).

For neutral conditions: Ch = 𝜅² / (log(Δz/z0b) * log(Δz/z0m)).
For stable/unstable: Ch = (u★ * φ_heat) / windspeed, where φ_heat is from universal functions.
"""
function heat_exchange_coefficient(
    param_set,
    L_MO,
    sc::Union{Fluxes, ValuesOnly, FluxesAndFrictionVelocity},
    scheme,
    tol_neutral,
    z0m,
    z0b,
)
    thermo_params = SFP.thermodynamics_params(param_set)
    𝜅 = SFP.von_karman_const(param_set)
    grav = SFP.grav(param_set)
    if abs(ΔDSEᵥ(param_set, sc)) <= tol_neutral
        Ch = 𝜅^2 / (log(Δz(sc) / z0b) * log(Δz(sc) / z0m))
    else
        ϕ_heat = compute_physical_scale_coeff(
            param_set,
            sc,
            L_MO,
            UF.HeatTransport(),
            scheme,
            z0m,
            z0b,
        )
        ustar = compute_ustar(param_set, L_MO, sc,  scheme, z0m, z0b)
        Ch = ustar * ϕ_heat / windspeed(sc)
    end
    return Ch
end

"""
    heat_exchange_coefficient(param_set, L_MO, sc,  scheme, tol_neutral, z0m, z0b)

Return Ch, the heat exchange coefficient given the
Monin-Obukhov lengthscale.
"""
function heat_exchange_coefficient(
    param_set,
    L_MO,
    sc::Coefficients,
    scheme,
    tol_neutral,
    z0m,
    z0b,
)
    return sc.Ch
end

"""
    momentum_fluxes(param_set, Cd, sc, scheme)

Compute momentum fluxes.

# Arguments
- `param_set`: Parameter set with physical constants
- `Cd`: Momentum exchange coefficient
- `sc`: Surface conditions container
- `scheme`: Discretization scheme

# Returns
Returns `(ρτxz, ρτyz)` tuple:
  - `ρτxz`: Eastward momentum flux (kg/m/s²)
  - `ρτyz`: Northward momentum flux (kg/m/s²)

Computed as: ρτ = -ρ_sfc * Cd * Δu * windspeed
"""
function momentum_fluxes(param_set, Cd, sc::AbstractSurfaceConditions, scheme)
    thermo_params = SFP.thermodynamics_params(param_set)
    ρ_sfc = TD.air_density(thermo_params, ts_sfc(sc))
    ρτxz = -ρ_sfc * Cd * Δu1(sc) * windspeed(sc)
    ρτyz = -ρ_sfc * Cd * Δu2(sc) * windspeed(sc)
    return (ρτxz, ρτyz)
end

"""
    sensible_heat_flux(param_set, Ch, sc, scheme)

Compute sensible heat flux.

# Arguments
- `param_set`: Parameter set with physical constants
- `Ch`: Heat exchange coefficient
- `sc`: Surface conditions container
- `scheme`: Discretization scheme

# Returns
Sensible heat flux (W/m²).

For `ValuesOnly` or `Coefficients`: Computed from dry static energy difference and evaporation.
"""
function sensible_heat_flux(
    param_set,
    Ch,
    sc::Union{ValuesOnly, Coefficients},
    scheme,
)
    thermo_params = SFP.thermodynamics_params(param_set)
    grav = SFP.grav(param_set)
    T_0 = SFP.T_0(param_set)
    LH_v0 = SFP.LH_v0(param_set)
    cp_m_in = TD.cp_m(thermo_params, ts_in(sc))
    cp_m_sfc = TD.cp_m(thermo_params, ts_sfc(sc))
    ρ_sfc = TD.air_density(thermo_params, ts_sfc(sc))
    T_in = TD.air_temperature(thermo_params, ts_in(sc))
    T_sfc = TD.air_temperature(thermo_params, ts_sfc(sc))
    hv_sfc = TD.specific_enthalpy_vapor(thermo_params, T_sfc)
    ΔΦ = grav * Δz(sc)
    ΔDSE = cp_m_in * (T_in - T_0) - cp_m_sfc * (T_sfc - T_0) + ΔΦ
    Φ_sfc = grav * z_sfc(sc)
    E = evaporation(param_set, sc, Ch)
    return -ρ_sfc * Ch * windspeed(sc) * ΔDSE + (hv_sfc + Φ_sfc - LH_v0) * E
end

"""
    sensible_heat_flux(param_set, Ch, sc, scheme)

Return known sensible heat flux from Fluxes or FluxesAndFrictionVelocity.
"""
function sensible_heat_flux(
    param_set,
    Ch,
    sc::Union{Fluxes, FluxesAndFrictionVelocity},
    scheme,
)
    return sc.shf
end

"""
    latent_heat_flux(param_set, L_MO, sc, scheme)

Return known latent heat flux from Fluxes or FluxesAndFrictionVelocity.
"""
function latent_heat_flux(
    param_set,
    L_MO,
    sc::Union{Fluxes, FluxesAndFrictionVelocity},
    scheme,
)
    return sc.lhf
end

"""
    latent_heat_flux(param_set, Ch, sc, scheme)

Compute latent heat flux.

# Arguments
- `param_set`: Parameter set with physical constants
- `Ch`: Heat exchange coefficient
- `sc`: Surface conditions container
- `scheme`: Discretization scheme

# Returns
Latent heat flux (W/m²).

Computed as: lhf = LH_v0 * evaporation(sc, Ch), where evaporation depends on specific humidity difference.
"""
function latent_heat_flux(
    param_set,
    Ch,
    sc::Union{ValuesOnly, Coefficients},
    scheme,
)
    LH_v0 = SFP.LH_v0(param_set)
    E = evaporation(param_set, sc, Ch)
    lhf = LH_v0 * E
    return lhf
end

"""
    evaporation(param_set, sc, Ch)

Compute evaporation rate.

# Arguments
- `param_set`: Parameter set with physical constants
- `sc`: Surface conditions container
- `Ch`: Heat exchange coefficient

# Returns
Evaporation rate (kg/m²/s).

For `Fluxes` or `FluxesAndFrictionVelocity`: E = lhf / LH_v0.
For `ValuesOnly` or `Coefficients`: E = -ρ_sfc * Ch * windspeed * Δqt * beta.
"""
function evaporation(
    param_set,
    sc::Union{Fluxes, FluxesAndFrictionVelocity},
    Ch,
)
    LH_v0 = SFP.LH_v0(param_set)
    return sc.lhf / LH_v0
end

"""
    evaporation(param_set, sc, Ch)

Compute evaporation rate using beta resistance factor.
"""
function evaporation(param_set, sc::Union{ValuesOnly, Coefficients}, Ch)
    thermo_params = SFP.thermodynamics_params(param_set)
    ρ_sfc = TD.air_density(thermo_params, ts_sfc(sc))
    return -ρ_sfc * Ch * windspeed(sc) * Δqt(param_set, sc) * sc.beta
end


"""
    compute_physical_scale_coeff(param_set, sc, L_MO, transport,  scheme, z0m, z0b)

Compute physical scale coefficient using finite volume scheme (Nishizawa 2018).

# Arguments
- `param_set`: Parameter set with physical constants
- `sc`: Surface conditions container
- `L_MO`: Monin-Obukhov length (m)
- `transport`: Transport type (MomentumTransport or HeatTransport)
- `uft`: Universal function type
- `scheme`: LayerAverageScheme (finite volume)
- `z0m`: Momentum roughness length (m)
- `z0b`: Scalar roughness length (m)

# Returns
Physical scale coefficient for momentum or heat transport (dimensionless).

Based on integrated universal functions for finite volume discretization.
"""
function compute_physical_scale_coeff(
    param_set::APS,
    sc::Union{ValuesOnly, Fluxes, FluxesAndFrictionVelocity},
    L_MO,
    transport,
    scheme::LayerAverageScheme,
    z0m,
    z0b,
)
    𝜅 = SFP.von_karman_const(param_set)
    π_group = UF.π_group(SFP.uf_params(param_set), transport)

    # Determine which z0 to use based on transport type
    z0_val = transport isa UF.MomentumTransport ? z0m : z0b

    R_z0 = 1 - z0_val / Δz(sc)
    denom1 = log(Δz(sc) / z0_val)
    denom2 = -UF.Psi(SFP.uf_params(param_set), Δz(sc) / L_MO, transport)
    denom3 = z0_val / Δz(sc) * UF.Psi(SFP.uf_params(param_set), z0_val / L_MO, transport)
    denom4 = R_z0 * (UF.psi(SFP.uf_params(param_set), z0_val / L_MO, transport) - 1)
    Σterms = denom1 + denom2 + denom3 + denom4
    return 𝜅 / (π_group * Σterms)
end


"""
    compute_physical_scale_coeff(param_set, sc, L_MO, transport,  scheme, z0m, z0b)

Compute physical scale coefficient using finite difference scheme (Byun 1990).

# Arguments
- `param_set`: Parameter set with physical constants
- `sc`: Surface conditions container
- `L_MO`: Monin-Obukhov length (m)
- `transport`: Transport type (MomentumTransport or HeatTransport)
- `uft`: Universal function type
- `scheme`: PointValueScheme (finite difference)
- `z0m`: Momentum roughness length (m)
- `z0b`: Scalar roughness length (m)

# Returns
Physical scale coefficient for momentum or heat transport (dimensionless).

Based on universal functions for finite difference discretization.
"""
function compute_physical_scale_coeff(
    param_set,
    sc::Union{ValuesOnly, Fluxes, FluxesAndFrictionVelocity},
    L_MO,
    transport,
    scheme::PointValueScheme,
    z0m,
    z0b,
)
    𝜅 = SFP.von_karman_const(param_set)
    π_group = UF.π_group(SFP.uf_params(param_set), transport)
    # Determine which z0 to use based on transport type
    z0_val = transport isa UF.MomentumTransport ? z0m : z0b
    denom1 = log(Δz(sc) / z0_val)
    denom2 = -UF.psi(SFP.uf_params(param_set), Δz(sc) / L_MO, transport)
    denom3 = UF.psi(SFP.uf_params(param_set), z0_val / L_MO, transport)
    Σterms = denom1 + denom2 + denom3
    return 𝜅 / (π_group * Σterms)
end

"""
    recover_profile(param_set, sc, L_MO, Z, X_star, X_sfc, transport, scheme)

Recover profiles using similarity theory (Nishizawa eq. 21, 22).

# Arguments
- `param_set`: Parameter set with physical constants
- `sc`: Surface conditions container
- `L_MO`: Monin-Obukhov length (m)
- `Z`: Height coordinate(s) within surface layer (m)
- `X_star`: Scale parameter for variable X
- `X_sfc`: Surface value for variable X
- `transport`: Transport type (MomentumTransport or HeatTransport)
- `scheme`: Discretization scheme

# Returns
Reconstructed profile values at heights Z.

# TODO: add tests
"""
function recover_profile(
    param_set::APS,
    sc::AbstractSurfaceConditions,
    L_MO,
    Z,
    X_star,
    X_sfc,
    transport,
    scheme::Union{LayerAverageScheme, PointValueScheme},
)
    𝜅 = SFP.von_karman_const(param_set)
    num1 = log(Z / z0(sc, transport))
    num2 = -UF.psi(SFP.uf_params(param_set), Z / L_MO, transport)
    num3 = UF.psi(SFP.uf_params(param_set), z0(sc, transport) / L_MO, transport)
    Σnum = num1 + num2 + num3
    return Σnum * X_star / 𝜅 + X_sfc
end

# For backwards compatibility with package extensions
if !isdefined(Base, :get_extension)
    include(joinpath("..", "ext", "CreateParametersExt.jl"))
end

end
