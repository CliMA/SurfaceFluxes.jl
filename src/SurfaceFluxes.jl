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

import KernelAbstractions
const KA = KernelAbstractions

using DocStringExtensions
const DSE = DocStringExtensions

import Thermodynamics
import Thermodynamics.ThermodynamicsParameters
const TD = Thermodynamics

import CLIMAParameters

import StaticArrays
const SA = StaticArrays

include("UniversalFunctions.jl")
import .UniversalFunctions
import .UniversalFunctions.AbstractUniversalFunctionParameters
const UF = UniversalFunctions


struct SurfaceFluxesParameters{FT, AUFPS <: AbstractUniversalFunctionParameters}
    von_karman_const::FT
    grav::FT
    gas_constant::FT
    kappa_d::FT
    T_0::FT
    molmass_dryair::FT
    molmass_water::FT
    cp_v::FT
    LH_v0::FT
    R_d::FT
    R_v::FT
    molmass_ratio::FT
    cp_d::FT
    cv_d::FT
    UFPS::AUFPS
    TPS::ThermodynamicsParameters
end

function SurfaceFluxesParameters(
    param_set,
    UFPS::AUFPS,
    TPS::ThermodynamicsParameters{FT}
) where {FT, AUFPS <: AbstractUniversalFunctionParameters}
    
    aliases = ["von_karman_const","grav","molmass_dryair",
        "molmass_water","gas_constant","kappa_d","T_0", "cp_v", "LH_v0"]
    (von_karman_const, grav, molmass_dryair, molmass_water,
    gas_constant, kappa_d, T_0, cp_v, LH_v0) = CLIMAParameters.get_parameter_values(
        param_set,
        aliases,
        "SurfaceFluxesParameters"
    )
    # derived parameters
    R_d = gas_constant / molmass_dryair
    R_v = gas_constant / molmass_water
    molmass_ratio = molmass_dryair / molmass_water
    cp_d = R_d / kappa_d
    cv_d = cp_d - R_d
    

    return SurfaceFluxesParameters{
        CLIMAParameters.get_parametric_type(param_set),
    }(
        von_karman_const,
        grav,
        gas_constant,
        kappa_d,
        T_0,
        molmass_dryair,
        molmass_water,
        cp_v,
        LH_v0,
        R_d,
        R_v,
        molmass_ratio,
        cp_d,
        cv_d,
        UFPS,
        TPS,
    )
end



    



abstract type SolverScheme end
struct FVScheme <: SolverScheme end
struct FDScheme <: SolverScheme end

# Allow users to skip error on non-convergence
# by importing:
# ```julia
# import SurfaceFluxes
# SurfaceFluxes.error_on_non_convergence() = false
# ```
# Error on convergence must be the default
# behavior because this can result in printing
# very large logs resulting in CI to seemingly hang.
error_on_non_convergence() = true

"""
    SurfaceFluxConditions

Surface flux conditions, returned from `surface_conditions`.

# Fields

$(DSE.FIELDS)
"""
struct SurfaceFluxConditions{FT <: AbstractFloat}
    L_MO::FT
    shf::FT
    lhf::FT
    buoy_flux::FT
    ρτxz::FT
    ρτyz::FT
    ustar::FT
    Cd::FT
    Ch::FT
end

function Base.show(io::IO, sfc::SurfaceFluxConditions)
    println(io, "----------------------- SurfaceFluxConditions")
    println(io, "L_MO                   = ", sfc.L_MO)
    println(io, "Sensible Heat Flux     = ", sfc.shf)
    println(io, "Latent Heat Flux       = ", sfc.lhf)
    println(io, "Buoyancy Flux          = ", sfc.buoy_flux)
    println(io, "Friction velocity u⋆   = ", sfc.ustar)
    println(io, "C_drag                 = ", sfc.Cd)
    println(io, "C_heat                 = ", sfc.Ch)
    println(io, "-----------------------")
end

"""
    SurfaceValues

Input container for state variables at the ground level.

# Fields

$(DSE.FIELDS)
"""
struct SurfaceValues{FT <: AbstractFloat, A, TS <: TD.ThermodynamicState}
    z::FT
    u::A
    ts::TS
end

"""
    InteriorValues

Input container for state variables at the first interior node.

# Fields

$(DSE.FIELDS)
"""
struct InteriorValues{FT <: AbstractFloat, A, TS <: TD.ThermodynamicState}
    z::FT
    u::A
    ts::TS
end

abstract type AbstractSurfaceConditions{FT <: AbstractFloat, VI <: InteriorValues, VS <: SurfaceValues} end


"""
    Fluxes

Input container for state variables, latent and sensible heat fluxes roughness lengths,
initial obukhov length and gustiness.

# Fields

$(DSE.FIELDS)
"""
Base.@kwdef struct Fluxes{FT, VI, VS} <: AbstractSurfaceConditions{FT, VI, VS}
    state_in::VI
    state_sfc::VS
    shf::FT
    lhf::FT
    z0m::FT
    z0b::FT
    L_MO_init::FT = FT(100)
    gustiness::FT = FT(1)
end

function Fluxes{FT}(; state_in, state_sfc, kwargs...) where {FT}
    types = (FT, typeof(state_in), typeof(state_sfc))
    return Fluxes{types...}(; state_in, state_sfc, kwargs...)
end


"""
    FluxesAndFrictionVelocity

Input container, given surface state variables, latent and sensible heat fluxes,
and the friction velocity, roughness lengths,
initial obukhov length and gustiness.

# Fields

$(DSE.FIELDS)
"""
Base.@kwdef struct FluxesAndFrictionVelocity{FT, VI, VS} <: AbstractSurfaceConditions{FT, VI, VS}
    state_in::VI
    state_sfc::VS
    shf::FT
    lhf::FT
    ustar::FT
    z0m::FT
    z0b::FT
    gustiness::FT = FT(1)
end

function FluxesAndFrictionVelocity{FT}(; state_in, state_sfc, kwargs...) where {FT}
    types = (FT, typeof(state_in), typeof(state_sfc))
    return FluxesAndFrictionVelocity{types...}(; state_in, state_sfc, kwargs...)
end

"""
    Coefficients

Input container, given surface state variables, and exchange coefficients,roughness lengths,
initial obukhov length and gustiness.

# Fields

$(DSE.FIELDS)
"""
Base.@kwdef struct Coefficients{FT, VI, VS} <: AbstractSurfaceConditions{FT, VI, VS}
    state_in::VI
    state_sfc::VS
    Cd::FT
    Ch::FT
    z0m::FT
    z0b::FT
    gustiness::FT = FT(1)
end

function Coefficients{FT}(; state_in, state_sfc, kwargs...) where {FT}
    types = (FT, typeof(state_in), typeof(state_sfc))
    return Coefficients{types...}(; state_in, state_sfc, kwargs...)
end


"""
    ValuesOnly

Input container, given only surface state variables, roughness lengths,
initial obukhov length and gustiness.

# Fields

$(DSE.FIELDS)
"""
Base.@kwdef struct ValuesOnly{FT, VI, VS} <: AbstractSurfaceConditions{FT, VI, VS}
    state_in::VI
    state_sfc::VS
    z0m::FT
    z0b::FT
    L_MO_init::FT = FT(100)
    gustiness::FT = FT(1)
end

function ValuesOnly{FT}(; state_in, state_sfc, kwargs...) where {FT}
    types = (FT, typeof(state_in), typeof(state_sfc))
    return ValuesOnly{types...}(; state_in, state_sfc, kwargs...)
end

ts_in(sc::AbstractSurfaceConditions) = sc.state_in.ts
ts_sfc(sc::AbstractSurfaceConditions) = sc.state_sfc.ts

z_in(sc::AbstractSurfaceConditions) = sc.state_in.z
z_sfc(sc::AbstractSurfaceConditions) = sc.state_sfc.z
Δz(sc::AbstractSurfaceConditions) = z_in(sc) - z_sfc(sc)

z0(sc::AbstractSurfaceConditions, ::UF.HeatTransport) = sc.z0b
z0(sc::AbstractSurfaceConditions, ::UF.MomentumTransport) = sc.z0m

Δu1(sc::AbstractSurfaceConditions) = sc.state_in.u[1] - sc.state_sfc.u[1]
Δu2(sc::AbstractSurfaceConditions) = sc.state_in.u[2] - sc.state_sfc.u[2]

qt_in(param_set::SurfaceFluxesParameters, sc::AbstractSurfaceConditions) = TD.total_specific_humidity(param_set.TPS, ts_in(sc))
qt_sfc(param_set::SurfaceFluxesParameters, sc::AbstractSurfaceConditions) = TD.total_specific_humidity(param_set.TPS, ts_sfc(sc))
Δqt(param_set::SurfaceFluxesParameters, sc::AbstractSurfaceConditions) = qt_in(param_set, sc) - qt_sfc(param_set, sc)

u_in(sc::AbstractSurfaceConditions) = sc.state_in.u
u_sfc(sc::AbstractSurfaceConditions) = sc.state_sfc.u

function windspeed(sc::AbstractSurfaceConditions)
    return max(hypot(Δu1(sc), Δu2(sc)), sc.gustiness)
end

"""
    surface_conditions(
        param_set::CLIMAParameters.AbstractParameterSet,
        sc::SurfaceFluxes.AbstractSurfaceConditions{FT},
        uft::UniversalFunctions.AbstractUniversalFunctionType,
        scheme::SurfaceFluxes.SolverScheme = FVScheme();
        tol::RootSolvers.AbstractTolerance =
            RootSolvers.SolutionTolerance(max(FT(1e-2), sqrt(eps(FT)))),
        maxiter::Int = 100,
    ) where {FT}

The main user facing function of the module.
It computes the surface conditions
based on the Monin-Obukhov similarity functions. Requires
information about thermodynamic parameters (`param_set`)
the surface state `sc`, the universal function type and
the discretisation `scheme`.

Result struct of type SurfaceFluxConditions{FT} contains:
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
    param_set::SurfaceFluxesParameters,
    sc::AbstractSurfaceConditions{FT},
    uft::UF.AUFT,
    scheme::SolverScheme = FVScheme();
    tol::RS.AbstractTolerance = RS.SolutionTolerance(max(FT(1e-2), sqrt(eps(FT)))),
    maxiter::Int = 100,
) where {FT}
    # Compute output variables for SurfaceFluxConditions
    L_MO = obukhov_length(param_set, sc, uft, scheme; tol, maxiter)
    ustar = compute_ustar(param_set, L_MO, sc, uft, scheme)
    Cd = momentum_exchange_coefficient(param_set, L_MO, sc, uft, scheme)
    Ch = heat_exchange_coefficient(param_set, L_MO, sc, uft, scheme)
    shf = sensible_heat_flux(param_set, Ch, sc, scheme)
    lhf = latent_heat_flux(param_set, Ch, sc, scheme)
    buoy_flux = compute_buoyancy_flux(param_set, shf, lhf, ts_in(sc), ts_sfc(sc), scheme)
    ρτxz, ρτyz = momentum_fluxes(param_set, Cd, sc, scheme)
    return SurfaceFluxConditions{FT}(L_MO, shf, lhf, buoy_flux, ρτxz, ρτyz, ustar, Cd, Ch)
end

"""
    obukhov_length(sfc::SurfaceFluxConditions)

    obukhov_length( # internal method
        param_set::AbstractParameterSet,
        sc::AbstractSurfaceConditions,
        uft,
        scheme;
        tol::RootSolvers.AbstractTolerance =
             RootSolvers.SolutionTolerance(max(1e-2, sqrt(eps))),
        maxiter::Int = 100
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
function obukhov_length end

obukhov_length(sfc::SurfaceFluxConditions) = sfc.L_MO

function obukhov_length(
    param_set,
    sc::AbstractSurfaceConditions{FT},
    uft::UF.AUFT,
    scheme;
    tol::RS.AbstractTolerance = RS.SolutionTolerance(max(FT(1e-2), sqrt(eps(FT)))),
    maxiter::Int = 100,
) where {FT}

    function root_l_mo(x_lmo)
        root = x_lmo - local_lmo(param_set, x_lmo, sc, uft, scheme)
        return root
    end
    sol = RS.find_zero(root_l_mo, RS.NewtonsMethodAD(sc.L_MO_init), RS.CompactSolution(), tol, maxiter)
    if sol.converged
        L_MO = sol.root
    else
        if error_on_non_convergence()
            KA.@print("-----------------------------------------\n")
            KA.@print("maxiter reached in SurfaceFluxes.jl:\n")
            KA.@print(", T_in = ", TD.air_temperature(param_set.TPS, ts_in(sc)))
            KA.@print(", T_sfc = ", TD.air_temperature(param_set.TPS, ts_sfc(sc)))
            KA.@print(", q_in = ", TD.total_specific_humidity(param_set.TPS, ts_in(sc)))
            KA.@print(", q_sfc = ", TD.total_specific_humidity(param_set.TPS, ts_sfc(sc)))
            KA.@print(", u_in = ", u_in(sc))
            KA.@print(", u_sfc = ", u_sfc(sc))
            error("Unconverged Surface Fluxes.")
        else
            # KA.@print("Warning: Unconverged Surface Fluxes. Returning neutral condition solution \n")
            L_MO = sol.root
        end
    end
    return L_MO
end

function obukhov_length(param_set::SurfaceFluxesParameters, sc::FluxesAndFrictionVelocity{FT}, uft::UF.AUFT, scheme; kwargs...) where {FT}
    return -sc.ustar^3 / FT(param_set.von_karman_const) / compute_buoyancy_flux(param_set, sc, scheme)
end

function obukhov_length(param_set::SurfaceFluxesParameters, sc::Coefficients{FT}, uft::UF.AUFT, scheme; kwargs...) where {FT}
    lhf = latent_heat_flux(param_set, sc.Ch, sc, scheme)
    shf = sensible_heat_flux(param_set, sc.Ch, sc, scheme)
    ustar = sqrt(sc.Cd) * windspeed(sc)
    buoyancy_flux = compute_buoyancy_flux(param_set, shf, lhf, ts_in(sc), ts_sfc(sc), scheme)
    return -ustar^3 / FT(param_set.von_karman_const) / buoyancy_flux
end

"""
    local_lmo(param_set, x_lmo, sc, uft, scheme)

Helper function for the iterative solver with the Monin-Obukhov length equation with buoyancy flux.
"""
function local_lmo(param_set::SurfaceFluxesParameters, x_lmo, sc::Fluxes, uft::UF.AUFT, scheme)
    κ = param_set.von_karman_const
    u_scale = compute_physical_scale_coeff(param_set, sc, x_lmo, UF.MomentumTransport(), uft, scheme)
    return -(windspeed(sc) * u_scale)^3 / κ / compute_buoyancy_flux(param_set, sc, scheme)
end

"""
    local_lmo(param_set, x_lmo, sc, uft, scheme)

Helper function for the iterative solver with the Monin-Obukhov
length equation with buoyancy star.
"""
function local_lmo(param_set::SurfaceFluxesParameters, x_lmo, sc::ValuesOnly, uft::UF.AUFT, scheme)
    κ = param_set.von_karman_const
    u_scale = compute_physical_scale_coeff(param_set, sc, x_lmo, UF.MomentumTransport(), uft, scheme)
    return (windspeed(sc) * u_scale)^2 / κ / compute_bstar(param_set, x_lmo, sc, uft, scheme)
end

"""
    compute_buoyancy_flux(param_set, shf, lhf, ts_in, ts_sfc, scheme)

Returns the buoyancy flux when the surface fluxes are known.
"""
function compute_buoyancy_flux(param_set::SurfaceFluxesParameters, shf::FT, lhf::FT, ts_in, ts_sfc, scheme) where {FT}

    grav = param_set.grav
    ε_vd = param_set.molmass_ratio

    cp_m = TD.cp_m(param_set.TPS, ts_in)
    L_v = TD.latent_heat_vapor(param_set.TPS, ts_in)
    ρ_sfc = TD.air_density(param_set.TPS, ts_sfc)
    T_in = TD.air_temperature(param_set.TPS, ts_in)
    return grav / ρ_sfc * (shf / cp_m / T_in + (ε_vd - 1) * lhf / L_v)
end

function compute_buoyancy_flux(param_set::SurfaceFluxesParameters, sc::Union{FluxesAndFrictionVelocity, Fluxes}, scheme)
    return compute_buoyancy_flux(param_set, sc.shf, sc.lhf, ts_in(sc), ts_sfc(sc), scheme)
end

"""
    compute_bstar(param_set, L_MO, sc, uft, scheme)

Returns buoyancy star based on known air densities.
"""
function compute_bstar(param_set::SurfaceFluxesParameters, L_MO, sc::AbstractSurfaceConditions{FT}, uft::UF.AUFT, scheme) where {FT}
    uf = UF.universal_func(uft, L_MO, param_set.UFPS)
    grav = param_set.grav

    ρ_in = TD.air_density(param_set.TPS, ts_in(sc))
    ρ_sfc = TD.air_density(param_set.TPS, ts_sfc(sc))

    ρ_star = compute_physical_scale_coeff(param_set, sc, L_MO, UF.HeatTransport(), uft, scheme) * (ρ_in - ρ_sfc)

    return -grav / ρ_in * ρ_star
end

"""
    compute_bstar(param_set, L_MO, sc, uft, scheme)

Returns buoyancy star based on known friction velocity  and fluxes.
"""
compute_bstar(param_set::SurfaceFluxesParameters, L_MO, sc::FluxesAndFrictionVelocity, uft, scheme) =
    -compute_buoyancy_flux(param_set, sc, scheme) / sc.ustar

"""
    compute_bstar(param_set, L_MO, sc, uft, scheme)

Return buoyancy star based on known fluxes.
"""
compute_bstar(param_set::SurfaceFluxesParameters, L_MO, sc::Fluxes, uft, scheme) =
    -compute_buoyancy_flux(param_set, sc, scheme) / compute_ustar(param_set, L_MO, sc, uft, scheme)


#=
    compute_ustar(
        param_set::AbstractParameterSet,
        L_MO,
        sc::AbstractSufaceCondition,
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
=#
function compute_ustar end

compute_ustar(param_set::SurfaceFluxesParameters, L_MO, sc::FluxesAndFrictionVelocity, uft, scheme) = sc.ustar

compute_ustar(param_set::SurfaceFluxesParameters, L_MO, sc::Fluxes, uft, scheme) =
    windspeed(sc) * compute_physical_scale_coeff(param_set, sc, L_MO, UF.MomentumTransport(), uft, scheme)

compute_ustar(param_set::SurfaceFluxesParameters, L_MO, sc::Coefficients, uft, scheme) = sqrt(sc.Cd) * (windspeed(sc))

compute_ustar(param_set::SurfaceFluxesParameters, L_MO, sc::ValuesOnly, uft, scheme) =
    windspeed(sc) * compute_physical_scale_coeff(param_set, sc, L_MO, UF.MomentumTransport(), uft, scheme)

"""
    momentum_exchange_coefficient(param_set, L_MO, sc, uft, scheme)

Compute and return Cd, the momentum exchange coefficient, given the
Monin-Obukhov lengthscale.
"""
function momentum_exchange_coefficient(
    param_set,
    L_MO,
    sc::Union{Fluxes, ValuesOnly, FluxesAndFrictionVelocity},
    uft::UF.AUFT,
    scheme,
)
    κ = param_set.von_karman_const
    transport = UF.MomentumTransport()
    FT = eltype(L_MO)
    ρ_in = TD.air_density(param_set.TPS, ts_in(sc))
    ρ_sfc = TD.air_density(param_set.TPS, ts_sfc(sc))
    Δρ = ρ_in - ρ_sfc
    ustar = compute_ustar(param_set, L_MO, sc, uft, scheme)
    if Δρ ≈ FT(0)
        Cd = (κ / log(Δz(sc) / z0(sc, transport)))^2
    else
        Cd = ustar^2 / windspeed(sc)^2
    end
    return Cd
end

"""
    momentum_exchange_coefficient(param_set, L_MO, sc, uft, scheme)

Return Cd, the momentum exchange coefficient.
"""
function momentum_exchange_coefficient(param_set::SurfaceFluxesParameters, L_MO, sc::Coefficients, uft, scheme)
    return sc.Cd
end

"""
    heat_exchange_coefficient(param_set, L_MO, sc, uft, scheme)

Compute and return Ch, the heat exchange coefficient given the
Monin-Obukhov lengthscale.
"""
function heat_exchange_coefficient(
    param_set,
    L_MO,
    sc::Union{Fluxes, ValuesOnly, FluxesAndFrictionVelocity},
    uft,
    scheme,
)
    FT = eltype(L_MO)
    Δθ = TD.dry_pottemp(param_set.TPS, ts_in(sc)) - TD.dry_pottemp(param_set.TPS, ts_sfc(sc))
    ustar = compute_ustar(param_set, L_MO, sc, uft, scheme)
    θstar = Δθ * compute_physical_scale_coeff(param_set, sc, L_MO, UF.HeatTransport(), uft, scheme)
    if Δθ ≈ FT(0) || windspeed(sc) ≈ FT(0)
        Ch = FT(0)
    else
        Ch = ustar * θstar / windspeed(sc) / Δθ
    end
    return Ch
end

"""
    heat_exchange_coefficient(param_set, L_MO, sc, uft, scheme)

Return Ch, the heat exchange coefficient given the
Monin-Obukhov lengthscale.
"""
function heat_exchange_coefficient(param_set::SurfaceFluxesParameters, L_MO, sc::Coefficients, uft, scheme)
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
function momentum_fluxes(param_set::SurfaceFluxesParameters, Cd, sc::AbstractSurfaceConditions, scheme)
    ρ_sfc = TD.air_density(param_set.TPS, ts_sfc(sc))
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
function sensible_heat_flux(param_set::SurfaceFluxesParameters, Ch::FT, sc::Union{ValuesOnly, Coefficients}, scheme) where {FT}
    grav = param_set.grav
    cp_d = param_set.cp_d
    R_d = param_set.R_d
    T_0 = param_set.T_0
    cp_m = TD.cp_m(param_set.TPS, ts_in(sc))
    ρ_sfc = TD.air_density(param_set.TPS, ts_sfc(sc))
    T_in = TD.air_temperature(param_set.TPS, ts_in(sc))
    T_sfc = TD.air_temperature(param_set.TPS, ts_sfc(sc))
    ΔT = T_in - T_sfc
    hd_sfc = cp_d * (T_sfc - T_0) + R_d * T_0
    ΔΦ = grav * Δz(sc)
    E = evaporation(sc, param_set, Ch)
    return -ρ_sfc * Ch * windspeed(sc) * (cp_m * ΔT + ΔΦ) - (hd_sfc) * E
end

function evaporation(sc::Union{ValuesOnly, Coefficients}, param_set, Ch)
    ρ_sfc = TD.air_density(param_set.TPS, ts_sfc(sc))
    return -ρ_sfc * Ch * windspeed(sc) * Δqt(param_set, sc)
end


"""
    sensible_heat_flux(param_set, Ch, sc, scheme)

In cases where surface fluxes are known,
return the known sensible heat flux.
"""
function sensible_heat_flux(param_set::SurfaceFluxesParameters, Ch, sc::Union{Fluxes, FluxesAndFrictionVelocity}, scheme)
    return sc.shf
end

"""
    latent_heat_flux(param_set, Ch, sc, scheme)

In cases where surface fluxes are known,
return the known latent heat flux.
"""
function latent_heat_flux(param_set::SurfaceFluxesParameters, L_MO, sc::Union{Fluxes, FluxesAndFrictionVelocity}, scheme)
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
function latent_heat_flux(param_set::SurfaceFluxesParameters, Ch::FT, sc::Union{ValuesOnly, Coefficients}, scheme) where {FT}
    grav::FT = param_set.grav
    ρ_sfc = TD.air_density(param_set.TPS, ts_sfc(sc))
    cp_v::FT = param_set.cp_v
    Lv_0::FT = param_set.LH_v0
    T_0::FT = param_set.T_0
    T_sfc = TD.air_temperature(param_set.TPS, ts_sfc(sc))
    hv_sfc = cp_v * (T_sfc - T_0) + Lv_0
    Φ_sfc = grav * z_sfc(sc)
    E = evaporation(sc, param_set, Ch)
    lhf = (hv_sfc + Φ_sfc) * E
    return lhf
end


"""
    compute_physical_scale_coeff(param_set, sc, L_MO, transport, uft, ::FVScheme)

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
    param_set::SurfaceFluxesParameters,
    sc::AbstractSurfaceConditions,
    L_MO::FT,
    transport,
    uft,
    ::FVScheme,
) where {FT}
    von_karman_const = param_set.von_karman_const
    uf = UF.universal_func(uft, L_MO,param_set.UFPS)
    _π_group = FT(UF.π_group(uf, transport))
    _π_group⁻¹ = (1 / _π_group)
    R_z0 = 1 - z0(sc, transport) / Δz(sc)
    denom1 = log(Δz(sc) / z0(sc, transport))
    denom2 = -UF.Psi(uf, Δz(sc) / uf.L, transport)
    denom3 = z0(sc, transport) / Δz(sc) * UF.Psi(uf, z0(sc, transport) / uf.L, transport)
    denom4 = R_z0 * (UF.psi(uf, z0(sc, transport) / uf.L, transport) - 1)
    Σterms = denom1 + denom2 + denom3 + denom4
    return _π_group⁻¹ * von_karman_const / Σterms
end


"""
    compute_physical_scale_coeff(param_set, sc, L_MO, transport, uft, ::FDScheme)

Computes the coefficient for the physical scale of a variable based on Byun(1990)
for the Finite Differneces scheme.

## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - L_MO: Monin-Obukhov length
  - transport: Transport type, (e.g. Momentum or Heat, used to determine physical scale coefficients)
  - uft: A Universal Function type, (returned by, e.g., Businger())
  - scheme: Discretization scheme (currently supports FD and FV)
"""
function compute_physical_scale_coeff(
    param_set,
    sc::AbstractSurfaceConditions,
    L_MO::FT,
    transport,
    uft,
    ::FDScheme,
) where {FT}
    von_karman_const = param_set.von_karman_const
    uf = UF.universal_func(uft, L_MO, param_set.UFPS)
    _π_group = FT(UF.π_group(uf, transport))
    _π_group⁻¹ = (1 / _π_group)
    denom1 = log(Δz(sc) / z0(sc, transport))
    denom2 = -UF.psi(uf, Δz(sc) / uf.L, transport)
    denom3 = UF.psi(uf, z0(sc, transport) / uf.L, transport)
    Σterms = denom1 + denom2 + denom3
    return _π_group⁻¹ * von_karman_const / Σterms
end

"""
    recover_profile(param_set, sc, L_MO, Z, X_in, X_sfc, transport, uft, scheme)

Recover profiles of variable X given values of Z coordinates. Follows Nishizawa equation (21,22)
## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - L_MO: Monin-Obukhov length
  - Z: Z coordinate(s) (within surface layer) for which variable values are required
  - X_in,X_sfc: For variable X, values at interior and surface nodes
  - transport: Transport type, (e.g. Momentum or Heat, used to determine physical scale coefficients)
  - uft: A Universal Function type, (returned by, e.g., Businger())
  - scheme: Discretization scheme (currently supports FD and FV)

# TODO: add tests
"""
function recover_profile(
    param_set::SurfaceFluxesParameters,
    sc::AbstractSurfaceConditions,
    L_MO,
    Z,
    X_in,
    X_sfc,
    transport,
    uft::UF.AUFT,
    scheme::Union{FVScheme, FDScheme},
) where {FT}
    @assert isless.(Z, sc.vals_in.z)
    uf = UF.universal_func(uft, L_MO, param_set.UFPS)
    von_karman_const::FT = param_set.von_karman_const
    _π_group = FT(UF.π_group(uf, transport))
    _π_group⁻¹ = (1 / _π_group)
    num1 = log(Z / z0(sc, transport))
    num2 = -UF.psi(uf, Z / L_MO, transport)
    num3 = UF.psi(uf, z0(sc, transport) / L_MO, transport)
    Σnum = num1 + num2 + num3
    ΔX = X_in - X_sfc
    return Σnum * compute_physical_scale_coeff(param_set, sc, L_MO, transport, uft, scheme) * _π_group⁻¹ * ΔX + X_sfc
end


end # SurfaceFluxes module
