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

include("UniversalFunctions.jl")
include("Parameters.jl")

import RootSolvers
const RS = RootSolvers

using DocStringExtensions
const DSE = DocStringExtensions

import Thermodynamics
const TD = Thermodynamics

import .UniversalFunctions
const UF = UniversalFunctions

import .Parameters

const SFP = Parameters
const APS = SFP.AbstractSurfaceFluxesParameters

abstract type SolverScheme end
struct LayerAverageScheme <: SolverScheme end
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

abstract type RoughnessModel end
struct CharnockRoughness <: RoughnessModel end
struct ScalarRoughness <: RoughnessModel end

abstract type AbstractSurfaceConditions{
    FT <: Real,
    SVA <: StateValues,
    SVB <: StateValues,
    RM <: RoughnessModel
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
        roughness_model
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
        roughness_model
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
    roughness_model::RM = ScalarRoughness()
) where {SVA, SVB, FT, RM}
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
        roughness_model
    )
end

function windspeed(sc::AbstractSurfaceConditions)
    return max(hypot(Δu1(sc), Δu2(sc)), sc.gustiness)
end

### Utilitity functions for calculations of differences between
### atmospheric state properties at the first interior node and

# Thermodynamic States
ts_in(sc::AbstractSurfaceConditions) = sc.state_in.ts
ts_sfc(sc::AbstractSurfaceConditions) = sc.state_sfc.ts

# Near-surface layer depth
z_in(sc::AbstractSurfaceConditions) = sc.state_in.z
z_sfc(sc::AbstractSurfaceConditions) = sc.state_sfc.z
Δz(sc::AbstractSurfaceConditions) = z_in(sc) - z_sfc(sc)

# Velocity
Δu1(sc::AbstractSurfaceConditions) = sc.state_in.u[1] - sc.state_sfc.u[1]
Δu2(sc::AbstractSurfaceConditions) = sc.state_in.u[2] - sc.state_sfc.u[2]

# Total Specific Humidity
qt_in(param_set::APS, sc::AbstractSurfaceConditions) =
    TD.total_specific_humidity(SFP.thermodynamics_params(param_set), ts_in(sc))
qt_sfc(param_set::APS, sc::AbstractSurfaceConditions) =
    TD.total_specific_humidity(SFP.thermodynamics_params(param_set), ts_sfc(sc))
Δqt(param_set::APS, sc::AbstractSurfaceConditions) =
    qt_in(param_set, sc) - qt_sfc(param_set, sc)

# Virtual Potential Temperature
θᵥ_in(param_set::APS, sc::AbstractSurfaceConditions) =
    TD.virtual_pottemp(SFP.thermodynamics_params(param_set), ts_in(sc))
θᵥ_sfc(param_set::APS, sc::AbstractSurfaceConditions) =
    TD.virtual_pottemp(SFP.thermodynamics_params(param_set), ts_sfc(sc))
Δθᵥ(param_set::APS, sc::AbstractSurfaceConditions) =
    θᵥ_in(param_set, sc) - θᵥ_sfc(param_set, sc)

# Virtual Dry Static Energy
DSEᵥ_in(param_set::APS, sc::AbstractSurfaceConditions) =
    TD.virtual_dry_static_energy(SFP.thermodynamics_params(param_set), 
                                 ts_in(sc), 
                                 SFP.grav(param_set)*z_in(sc))
DSEᵥ_sfc(param_set::APS, sc::AbstractSurfaceConditions) =
    TD.virtual_dry_static_energy(SFP.thermodynamics_params(param_set), 
                                 ts_sfc(sc), 
                                 SFP.grav(param_set)*z_sfc(sc))
ΔDSEᵥ(param_set::APS, sc::AbstractSurfaceConditions) =
    DSEᵥ_in(param_set, sc) - DSEᵥ_sfc(param_set, sc)

# Roughness Model Computations
function compute_z0(u★, sfc_param_set, 
                    sc::AbstractSurfaceConditions, ::ScalarRoughness, ::UF.MomentumTransport)
    return sc.z0m
end

function compute_z0(u★, sfc_param_set, 
                    sc, ::CharnockRoughness, ::UF.MomentumTransport)
    𝛼 = eltype(u★)(0.011)
    𝑔 = SFP.grav(sfc_param_set)
    return 𝛼 * u★^2 / 𝑔
end
function compute_z0(u★, sfc_param_set, 
                    sc, ::CharnockRoughness, ::UF.HeatTransport)
    𝛼 = eltype(u★)(0.011)
    𝑔 = SFP.grav(sfc_param_set)
    return 𝛼 * u★^2 / 𝑔
end

function compute_z0(u★, sfc_param_set, 
                    sc::AbstractSurfaceConditions, ::ScalarRoughness, ::UF.HeatTransport)
    return sc.z0m
end
function compute_z0(u★, sfc_param_set, 
                    sc::Coefficients, ::ScalarRoughness, ::UF.MomentumTransport)
    return compute_z0(u★, sfc_param_set, sc, CharnockRoughness(), UF.MomentumTransport())
end
function compute_z0(u★, sfc_param_set, 
                    sc::Coefficients, ::ScalarRoughness, ::UF.HeatTransport)
    return compute_z0(u★, sfc_param_set, sc, CharnockRoughness(), UF.MomentumTransport())
end
function compute_z0(u★, sfc_param_set, 
                    sc, ::ScalarRoughness, ::UF.HeatTransport)
    return sc.z0b
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
        noniterative_stable_sol::Bool=true,
    )

The main user facing function of the module.
It computes the surface conditions
based on the Monin-Obukhov similarity functions. Requires
information about thermodynamic parameters (`param_set`)
the surface state `sc`, the universal function type and
the discretisation `scheme`. Default tolerance for
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
"""
function surface_conditions(
    param_set::APS{FT},
    sc::AbstractSurfaceConditions,
    scheme::SolverScheme = PointValueScheme();
    tol_neutral = SFP.cp_d(param_set) / 100,
    tol::RS.AbstractTolerance = RS.RelativeSolutionTolerance(FT(0.01)),
    maxiter::Int = 30,
    soltype::RS.SolutionType = RS.CompactSolution(),
    noniterative_stable_sol::Bool = true,
) where {FT}
    uft = SFP.universal_func_type(param_set)
    X★ = obukhov_similarity_solution(
        param_set,
        sc,
        uft,
        scheme,
        tol,
        tol_neutral,
        maxiter,
        soltype,
        noniterative_stable_sol,
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
        uft,
        scheme,
        tol_neutral,
    )
    Ch =
        heat_exchange_coefficient(
        param_set, 
        L_MO, 
        ustar, 
        sc, 
        uft, 
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
    tol::RS.AbstractTolerance = RS.RelativeSolutionTolerance(FT(0.01)),
    maxiter::Int = 30,
    soltype::RS.SolutionType = RS.CompactSolution(),
    noniterative_stable_sol::Bool = true,
) where {FT}
    uft = SFP.universal_func_type(param_set)
    L_MO, _ = obukhov_similarity_solution(
        param_set,
        sc,
        uft,
        scheme,
        tol,
        tol_neutral,
        maxiter,
        soltype,
        noniterative_stable_sol,
    )
    ustar = sc.ustar
    Cd = momentum_exchange_coefficient(param_set, L_MO, ustar, sc, uft, scheme, tol_neutral)
    Ch = heat_exchange_coefficient(param_set, L_MO, ustar, sc, uft, scheme, tol_neutral)
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
    sc::Coefficients,
    scheme::SolverScheme = PointValueScheme();
    tol_neutral = SFP.cp_d(param_set) / 100,
    tol::RS.AbstractTolerance = RS.RelativeSolutionTolerance(FT(0.01)),
    maxiter::Int = 30,
    soltype::RS.SolutionType = RS.CompactSolution(),
    noniterative_stable_sol::Bool = true,
) where {FT}
    uft = SFP.universal_func_type(param_set)
    ustar = sqrt(sc.Cd) * windspeed(sc)
    L_MO, _ = obukhov_similarity_solution(param_set, sc, uft, scheme)
    Cd = momentum_exchange_coefficient(
        param_set,
        nothing, 
        nothing,
        sc,
        uft,
        scheme,
        tol_neutral,
    )
    Ch =
        heat_exchange_coefficient(
        param_set, 
        nothing, 
        nothing, 
        sc, 
        uft, 
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
        soltype,
    )

Compute and return the Monin-Obukhov lengthscale (LMO).

The internal method for computing LMO depends on the
particular surface condition `sc`, of which there are
several options:

 - `FluxesAndFrictionVelocity`
 - `Coefficients`

## `AbstractSurfaceConditions` (fallback)

The Monin-Obukhov length is computed by solving a non-linear
equation given a tolerance `tol` and maximum iterations `maxiter`.

## `FluxesAndFrictionVelocity`

Surface fluxes and friction velocity are known.
Iterations are not needed to determine LMO.

## `Coefficients`

Exchange coefficients are known.
Iterations are not needed to determine LMO.
"""
function obukhov_similarity_solution end

obukhov_similarity_solution(sfc::SurfaceFluxConditions) = sfc.L_MO

function non_zero(v::FT) where {FT}
    sign_of_v = v == 0 ? 1 : sign(v)
    return abs(v) < eps(FT) ? eps(FT) * sign_of_v : v
end

function compute_Fₘₕ(sc, ufₛ, ζ, 𝓁, transport)
    ψ = UF.psi(ufₛ, ζ, transport)
    ψ₀ = UF.psi(ufₛ, 𝓁 * ζ / Δz(sc), transport)
    return log(Δz(sc)/𝓁) - ψ + ψ₀
end

function obukhov_similarity_solution(
    param_set::APS{FT},
    sc::Union{Fluxes, ValuesOnly},
    uft::UF.AUFT,
    scheme,
    tol,
    tol_neutral,
    maxiter,
    soltype,
    noniterative_stable_sol,
) where {FT}
    thermo_params = SFP.thermodynamics_params(param_set)
    ufparams = SFP.uf_params(param_set)
    grav = SFP.grav(param_set)
    δ = sign(ΔDSEᵥ(param_set, sc))
    if ΔDSEᵥ(param_set, sc) >= FT(0)
        X★₀ = SimilarityScaleVars(FT(δ), FT(δ), FT(δ), 
                                  FT(10),
                                  FT(0.0001), FT(0.0001), FT(0.0001))
        X★ = obukhov_iteration(X★₀, sc, uft, scheme, param_set)
        return  X★
    else
        X★₀ = SimilarityScaleVars(FT(δ), FT(δ), FT(δ), 
                                  FT(-10), 
                                  FT(0.0001), FT(0.0001), FT(0.0001))
        X★ = obukhov_iteration(X★₀,sc, uft, scheme, param_set)
        return  X★
    end
end

function obukhov_similarity_solution(
    param_set,
    sc::FluxesAndFrictionVelocity,
    uft::UF.AUFT,
    scheme,
    args...,
)
    return (-sc.ustar^3 / SFP.von_karman_const(param_set) /
            non_zero(compute_buoyancy_flux(param_set, sc, scheme)), sc.ustar)
end

function obukhov_similarity_solution(
    param_set,
    sc::Coefficients,
    uft::UF.AUFT,
    scheme,
    args...,
)
    lhf = latent_heat_flux(param_set, sc.Ch, sc, scheme)
    shf = sensible_heat_flux(param_set, sc.Ch, sc, scheme)
    ustar = sqrt(sc.Cd) * windspeed(sc)
    buoyancy_flux = compute_buoyancy_flux(
        param_set,
        shf,
        lhf,
        ts_in(sc),
        ts_sfc(sc),
        scheme,
    )
    return (-ustar^3 / SFP.von_karman_const(param_set) / non_zero(buoyancy_flux), ustar)
end

"""
    compute_buoyancy_flux(param_set, shf, lhf, ts_in, ts_sfc, scheme)

Returns the buoyancy flux when the surface fluxes are known.
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
        sc::AbstractSurfaceCondition,
        uft,
        scheme
    )

Return the friction velocity. This method is dispatched
by the surface condition:

## `sc::FluxesAndFrictionVelocity`

Friction velocity is known.

## `sc::Fluxes`

Compute given the Monin-Obukhov lengthscale.

## `sc::Coefficients`

Compute given the exchange coefficients.

## `sc::ValuesOnly`
Compute given the Monin-Obukhov lengthscale.
"""
function compute_ustar end

compute_ustar(param_set, L_MO, 𝓁, sc::FluxesAndFrictionVelocity, uft, scheme) =
    sc.ustar

compute_ustar(param_set, L_MO, 𝓁, sc::Fluxes, uft, scheme) =
    windspeed(sc) * compute_physical_scale_coeff(
        param_set,
        sc,
        L_MO,
        𝓁,
        UF.MomentumTransport(),
        uft,
        scheme,
    )

compute_ustar(param_set, L_MO, 𝓁, sc::Coefficients, uft, scheme) =
    sqrt(sc.Cd) * (windspeed(sc))

compute_ustar(param_set, L_MO, 𝓁, sc::ValuesOnly, uft, scheme) =
    windspeed(sc) * compute_physical_scale_coeff(
        param_set,
        sc,
        L_MO,
        𝓁,
        UF.MomentumTransport(),
        uft,
        scheme,
    )

"""
    momentum_exchange_coefficient(param_set, L_MO, sc, uft, scheme)

Compute and return Cd, the momentum exchange coefficient, given the
Monin-Obukhov lengthscale.
"""
function momentum_exchange_coefficient(
    param_set,
    L_MO,
    u★,
    sc::Union{Fluxes, ValuesOnly, FluxesAndFrictionVelocity},
    uft::UF.AUFT,
    scheme,
    tol_neutral,
)
    thermo_params = SFP.thermodynamics_params(param_set)
    κ = SFP.von_karman_const(param_set)
    grav = SFP.grav(param_set)
    transport = UF.MomentumTransport()
    𝓁 = compute_z0(u★, param_set, sc, sc.roughness_model, UF.MomentumTransport())
    if abs(ΔDSEᵥ(param_set, sc)) <= tol_neutral
        Cd = (κ / log(Δz(sc) / 𝓁))^2
    else
        ustar = compute_ustar(param_set, L_MO, 𝓁, sc, uft, scheme)
        Cd = ustar^2 / windspeed(sc)^2
    end
    return Cd
end

"""
    momentum_exchange_coefficient(param_set, L_MO, sc, uft, scheme, tol_neutral)

Return Cd, the momentum exchange coefficient.
"""
function momentum_exchange_coefficient(
    param_set,
    L_MO,
    u★,  
    sc::Coefficients,
    uft,
    scheme,
    tol_neutral,
)
    return sc.Cd
end

"""
    heat_exchange_coefficient(param_set, L_MO, sc, uft, scheme, tol_neutral)

Compute and return Ch, the heat exchange coefficient given the
Monin-Obukhov lengthscale.
"""
function heat_exchange_coefficient(
    param_set,
    L_MO,
    u★, 
    sc::Union{Fluxes, ValuesOnly, FluxesAndFrictionVelocity},
    uft,
    scheme,
    tol_neutral,
)       
    thermo_params = SFP.thermodynamics_params(param_set)
    transport = UF.HeatTransport()
    κ = SFP.von_karman_const(param_set)
    grav = SFP.grav(param_set)
    𝓁u = compute_z0(u★, param_set, sc, sc.roughness_model, UF.MomentumTransport())
    𝓁θ = compute_z0(u★, param_set, sc, sc.roughness_model, UF.HeatTransport())
    if abs(ΔDSEᵥ(param_set, sc)) <= tol_neutral
        Ch = κ^2 / (log(Δz(sc) / 𝓁θ) * log(Δz(sc) / 𝓁u))
    else
        ϕ_heat = compute_physical_scale_coeff(
            param_set,
            sc,
            L_MO, 
            𝓁θ,
            transport,
            uft,
            scheme,
        )
        ustar = compute_ustar(param_set, L_MO, 𝓁u, sc, uft, scheme)
        Ch = ustar * ϕ_heat / windspeed(sc)
    end
    return Ch
end

"""
    heat_exchange_coefficient(param_set, L_MO, sc, uft, scheme)

Return Ch, the heat exchange coefficient given the
Monin-Obukhov lengthscale.
"""
function heat_exchange_coefficient(
    param_set,
    L_MO,
    u★,
    sc::Coefficients,
    uft,
    scheme,
    tol_neutral,
)
    return sc.Ch
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
    sensible_heat_flux(param_set, Ch, sc, scheme)

Compute and return the sensible heat fluxes
## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - Ch: Thermal exchange coefficient
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - scheme: Discretization scheme (currently supports FD and FV)
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

In cases where surface fluxes are known,
return the known sensible heat flux.
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
    latent_heat_flux(param_set, Ch, sc, scheme)

In cases where surface fluxes are known,
return the known latent heat flux.
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

Compute and return the latent heat flux
## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - Ch: Thermal exchange coefficient
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - scheme: Discretization scheme (currently supports FD and FV)
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

Compute and return the evaporation. When `sc` is `Fluxes` or `FluxesAndFrictionVelocity`,
evaporation is directly calculated from the latent heat flux.
## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - Ch: Thermal exchange coefficient
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

Compute and return the evaporation. When `sc` is `ValuesOnly` or `Coefficients`, a `beta` factor
is used to represent the resistance of the surface.
## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - Ch: Thermal exchange coefficient
"""
function evaporation(param_set, sc::Union{ValuesOnly, Coefficients}, Ch)
    thermo_params = SFP.thermodynamics_params(param_set)
    ρ_sfc = TD.air_density(thermo_params, ts_sfc(sc))
    return -ρ_sfc * Ch * windspeed(sc) * Δqt(param_set, sc) * sc.beta
end


"""
    compute_physical_scale_coeff(param_set, sc, L_MO, transport, uft, ::LayerAverageScheme)

Computes the coefficient for the physical scale of a variable based on Nishizawa(2018)
for the FV scheme.

## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - L_MO: Monin-Obukhov length
  - transport: Transport type, (e.g. Momentum or Heat, used to determine physical scale coefficients)
  - uft: A Universal Function type, (returned by, e.g. Businger())
  - scheme: Discretization scheme (currently supports FD and FV)
"""
function compute_physical_scale_coeff(
    param_set::APS,
    sc::Union{ValuesOnly, Fluxes, FluxesAndFrictionVelocity},
    L_MO,
    𝓁, 
    transport,
    uft,
    ::LayerAverageScheme,
)
    von_karman_const = SFP.von_karman_const(param_set)
    uf = UF.universal_func(uft, L_MO, SFP.uf_params(param_set))
    π_group = UF.π_group(uf, transport)
    R_z0 = 1 - 𝓁 / Δz(sc)
    denom1 = log(Δz(sc) / 𝓁)
    denom2 = -UF.Psi(uf, Δz(sc) / uf.L, transport)
    denom3 =
        𝓁 / Δz(sc) *
        UF.Psi(uf, 𝓁 / uf.L, transport)
    denom4 = R_z0 * (UF.psi(uf, 𝓁 / uf.L, transport) - 1)
    Σterms = denom1 + denom2 + denom3 + denom4
    return von_karman_const / (π_group * Σterms)
end


"""
    compute_physical_scale_coeff(param_set, sc, L_MO, transport, uft, ::PointValueScheme)

Computes the coefficient for the physical scale of a variable based on Byun (1990)
for the Finite Differences scheme.

## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - L_MO: Monin-Obukhov length
  - transport: Transport type, (e.g. Momentum or Heat, used to determine physical scale coefficients)
  - uft: A Universal Function type, (returned by, e.g. Businger())
  - scheme: Discretization scheme (currently supports FD and FV)
"""
function compute_physical_scale_coeff(
    param_set,
    sc::Union{ValuesOnly, Fluxes, FluxesAndFrictionVelocity},
    L_MO,
    𝓁, 
    transport,
    uft,
    ::PointValueScheme,
)
    von_karman_const = SFP.von_karman_const(param_set)
    uf = UF.universal_func(uft, L_MO, SFP.uf_params(param_set))
    π_group = UF.π_group(uf, transport)
    denom1 = log(Δz(sc) / 𝓁)
    denom2 = -UF.psi(uf, Δz(sc) / uf.L, transport)
    denom3 = UF.psi(uf, 𝓁 / uf.L, transport)
    Σterms = denom1 + denom2 + denom3
    return von_karman_const / (π_group * Σterms)
end

"""
    recover_profile(param_set, sc, L_MO, Z, X_in, X_sfc, transport, scheme)

Recover profiles of variable X given values of Z coordinates. Follows Nishizawa equation (21,22)
## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - L_MO: Monin-Obukhov length
  - Z: Z coordinate(s) (within surface layer) for which variable values are required
  - X_star: Scale parameter for variable X
  - X_sfc: For variable X, values at surface nodes
  - transport: Transport type, (e.g. Momentum or Heat, used to determine physical scale coefficients)
  - scheme: Discretization scheme (currently supports FD and FV)

# TODO: add tests
"""
function recover_profile(
    param_set::APS,
    sc::AbstractSurfaceConditions,
    L_MO,
    𝓁,
    Z,
    X_star,
    X_sfc,
    transport,
    scheme::Union{LayerAverageScheme, PointValueScheme},
)
    ufp = SFP.uf_params(param_set)
    uft = UF.universal_func_type(typeof(ufp))
    uf = UF.universal_func(uft, L_MO, ufp)
    von_karman_const = SFP.von_karman_const(param_set)
    num1 = log(Z / 𝓁)
    num2 = -UF.psi(uf, Z / L_MO, transport)
    num3 = UF.psi(uf, 𝓁 / L_MO, transport)
    Σnum = num1 + num2 + num3
    return Σnum * X_star / von_karman_const + X_sfc
end

sfc_param_set(FT, UFT) = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
thermo_params(param_set) = SFP.thermodynamics_params(param_set)
uft(param_set) = UF.universal_func_type(param_set)

struct SimilarityScaleVars 
    u★
    DSEᵥ★
    q★
    L★
    𝓁u
    𝓁θ
    𝓁q
end

@inline function buoyancy_scale(DSEᵥ★, q★, thermo_params, 𝒬, 𝑔)
    𝒯ₐ = TD.virtual_temperature(thermo_params, 𝒬)
    qₐ = TD.vapor_specific_humidity(thermo_params, 𝒬)
    ε  = TD.Parameters.Rv_over_Rd(thermo_params)
    δ  = ε - 1 
    cp_v = TD.Parameters.cp_v(thermo_params)
    # Convert DSEᵥ★ (energy scale) to temperature scale for buoyancy calculation
    θ★_equiv = DSEᵥ★ / cp_v
    b★ = 𝑔 / 𝒯ₐ * (θ★_equiv * (1 + δ * qₐ) + δ * 𝒯ₐ * q★)
    return b★
end

"""
    iterate_interface_fluxes()
"""
function iterate_interface_fluxes(sc::Union{ValuesOnly, Fluxes},
                                  DSEᵥ₀, qₛ,
                                  approximate_interface_state,
                                  atmosphere_state,
                                  surface_state,
                                  uft::UF.AUFT,
                                  scheme::SolverScheme,
                                  param_set::APS)

    # Stability function type and problem parameters
    uft = SFP.universal_func_type(param_set)
    𝜅 = SFP.von_karman_const(param_set)
    ufparams = SFP.uf_params(param_set)
    𝑔 = SFP.grav(param_set)

    thermo_params = SFP.thermodynamics_params(param_set)

    #DSEᵥ₀ is passed in and represents surface DSEᵥ
    qₛ = qt_sfc(param_set, sc)
    ΔDSEᵥ_diff = ΔDSEᵥ(param_set, sc)
    Δq = Δqt(param_set, sc)

    # "Initial" approximate scales because we will recompute them
    # Get these from the model properties directly, given
    # some initial guess for L★ and u★
    u★ = approximate_interface_state.u★
    DSEᵥ★ = approximate_interface_state.DSEᵥ★
    q★ = Δq == eltype(𝑔)(0) ? approximate_interface_state.q★ : eltype(𝑔)(0)
    L★ = approximate_interface_state.L★
    𝓁u = compute_z0(u★, param_set, sc, sc.roughness_model, UF.MomentumTransport())
    𝓁θ = compute_z0(u★, param_set, sc, sc.roughness_model, UF.HeatTransport())
    𝓁q = compute_z0(u★, param_set, sc, sc.roughness_model, UF.HeatTransport())

    ζ = Δz(sc) / L★ # Stability Parameter

    # Stability functions for momentum, heat, and vapor
    ufp = SFP.uf_params(param_set)
    uf = UF.universal_func(uft, L★, ufp)

    # Compute interior and surface thermodynamic state
    𝒬ₛ = ts_sfc(sc)

    # Compute Monin--Obukhov length scale depending on a `buoyancy flux`
    b★ = buoyancy_scale(DSEᵥ★, q★, thermo_params, 𝒬ₛ, 𝑔)
    Uᴳ = sc.gustiness
    U = sqrt(windspeed(sc)^2)

    # Transfer coefficients at height `h`
    𝜅 = SFP.von_karman_const(param_set)
    δdseᵥ = ΔDSEᵥ(param_set, sc)
    L★ = ifelse(b★ == 0, sign(δdseᵥ) * Inf, u★^2 / (𝜅 * b★))
    ζ = Δz(sc) / L★      

    Pr = UF.Pr_0(uf)
    
    # Compute new values for the scale parameters given the relation
    χu = 𝜅 / compute_Fₘₕ(sc, uf, ζ, 𝓁u, UF.MomentumTransport())
    χθ = 𝜅 / Pr / compute_Fₘₕ(sc, uf, ζ, 𝓁θ, UF.HeatTransport())
    χq = 𝜅 / Pr / compute_Fₘₕ(sc, uf, ζ, 𝓁q, UF.HeatTransport())

    # Recompute
    u★ = χu * U
    DSEᵥ★ = χθ * ΔDSEᵥ_diff
    q★ = χq * Δq

    return SimilarityScaleVars(u★, DSEᵥ★, q★, L★, 𝓁u, 𝓁θ, 𝓁q)
end

function obukhov_iteration(X★, sc,
                           uft, scheme, param_set,
                           tol = sqrt(eps(eltype(X★.u★))), maxiter=30)
       DSEᵥ₀ = DSEᵥ_sfc(param_set, sc)
       q₀ = qt_sfc(param_set, sc)
       ts₀ = ts_sfc(sc)
       ts₁ = ts_in(sc)
       for ii = 1:maxiter
              X★₀ = X★
              X★ = iterate_interface_fluxes(sc,
                                          DSEᵥ₀, q₀,
                                          X★,
                                          ts₁,
                                          ts₀,
                                          uft,
                                          scheme,
                                          param_set)
              # Generalize and define a method for the norm evaluation
              if abs(X★.u★ - X★₀.u★) <= tol &&  
                 abs(X★.L★ - X★₀.L★) <= 100tol &&
                 abs(X★.DSEᵥ★ - X★₀.DSEᵥ★) <= tol 
                     return X★
              else
                     continue
              end
       end
       return X★
end

# For backwards compatibility with package extensions
if !isdefined(Base, :get_extension)
    include(joinpath("..", "ext", "CreateParametersExt.jl"))
end

end # SurfaceFluxes module
