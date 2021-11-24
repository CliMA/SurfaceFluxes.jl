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

import NonlinearSolvers
import RootSolvers
const NS = NonlinearSolvers
const RS = RootSolvers

import KernelAbstractions
const KA = KernelAbstractions

using DocStringExtensions
const DSE = DocStringExtensions

using Thermodynamics
const TD = Thermodynamics

import CLIMAParameters
const CP = CLIMAParameters
const CPP = CP.Planet
const APS = CP.AbstractEarthParameterSet
const CPSGS = CP.SubgridScale

import StaticArrays
const SA = StaticArrays

include("UniversalFunctions.jl")
import .UniversalFunctions
const UF = UniversalFunctions

abstract type SolverScheme end
struct FVScheme <: SolverScheme end
struct DGScheme <: SolverScheme end

export surface_conditions, exchange_coefficients, recover_profile, Values

"""
    SurfaceFluxConditions{FT}

Surface flux conditions, returned from `surface_conditions`.

# Fields

$(DSE.FIELDS)
"""
struct SurfaceFluxConditions{FT}
    L_MO::FT
    shf::FT
    lhf::FT
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
    println(io, "Friction velocity u⋆   = ", sfc.ustar)
    println(io, "C_drag                 = ", sfc.Cd)
    println(io, "C_heat                 = ", sfc.Ch)
    println(io, "-----------------------")
end

abstract type AbstractSurfaceLocation end

abstract type AbstractSurfaceConditions end

"""
    Fluxes{VI, VS, AUF, FT}

Input container, given surface state variables, latent and sensible heat fluxes.

# Fields

$(DSE.FIELDS)
"""
struct Fluxes{VI, VS, AUF, FT} <: AbstractSurfaceConditions
    state_in::VI
    state_sfc::VS
    shf::FT
    lhf::FT
    uf::AUF
    z0m::FT
    z0b::FT
    L_MO_init::FT
    gustiness::FT
end
function GivenFluxes(
    state_in::VI,
    state_sfc::VS,
    shf::FT,
    lhf::FT,
    uf::AUF,
    z0m::FT,
    z0b::FT;
    L_MO_init::FT = FT(100),
    gustiness::FT = FT(1),
) where {VI, VS, AUF, FT <: AbstractFloat}
    return Fluxes(
        state_in,
        state_sfc,
        shf,
        lhf,
        uf,
        z0m,
        z0b,
        L_MO_init,
        gustiness,
    )
end


"""
    FluxesAndFrictionVelocity{VI, VS, AUF, FT}

Input container, given surface state variables, latent and sensible heat fluxes,
and the friction velocity.

# Fields

$(DSE.FIELDS)
"""
struct FluxesAndFrictionVelocity{VI, VS, AUF, FT} <: AbstractSurfaceConditions
    state_in::VI
    state_sfc::VS
    shf::FT
    lhf::FT
    ustar::FT
    uf::AUF
    z0m::FT
    z0b::FT
    gustiness::FT
end
function GivenFluxesAndFrictionVelocity(
    state_in::VI,
    state_sfc::VS,
    shf::FT,
    lhf::FT,
    ustar::FT,
    uf::AUF,
    z0m::FT,
    z0b::FT;
    gustiness::FT = FT(1.0),
) where {VI, VS, AUF, FT <: AbstractFloat}
    return FluxesAndFrictionVelocity(
        state_in,
        state_sfc,
        shf,
        lhf,
        ustar,
        uf,
        z0m,
        z0b,
        gustiness,
    )
end

"""
    Coefficients{VI, VS, AUF, FT}

Input container, given surface state variables, and exchange coefficients.

# Fields

$(DSE.FIELDS)
"""
Base.@kwdef struct Coefficients{VI, VS, AUF, FT} <: AbstractSurfaceConditions
    state_in::VI
    state_sfc::VS
    Cd::FT
    Ch::FT
    uf::AUF
    z0m::FT
    z0b::FT
    gustiness::FT
end
function GivenCoefficients(
    state_in::VI,
    state_sfc::VS,
    Cd::FT,
    Ch::FT,
    uf::AUF,
    z0m::FT,
    z0b::FT;
    gustiness::FT = FT(1.0),
) where {VI, VS, AUF, FT <: AbstractFloat}
    return Coefficients(state_in, state_sfc, Cd, Ch, uf, z0m, z0b, gustiness)
end

"""
    ValuesOnly{VI, VS, AUF, FT}

Input container, given only surface state variables.

# Fields

$(DSE.FIELDS)
"""
struct ValuesOnly{VI, VS, AUF, FT} <: AbstractSurfaceConditions
    state_in::VI
    state_sfc::VS
    uf::AUF
    z0m::FT
    z0b::FT
    L_MO_init::FT
    gustiness::FT
end
function GivenValuesOnly(
    state_in::VI,
    state_sfc::VS,
    uf::AUF,
    z0m::FT,
    z0b::FT,
    L_MO_init::FT;
    gustiness::FT = FT(1.0),
) where {VI, VS, AUF, FT <: AbstractFloat}
    return ValuesOnly(state_in, state_sfc, uf, z0m, z0b, L_MO_init, gustiness)
end

"""
  ValuesSurface{FT, A, TS <: TD.ThermodynamicState}

Input container for state variables at the ground level. 

# Fields

$(DSE.FIELDS)
"""
struct ValuesSurface{FT, A, TS <: TD.ThermodynamicState} <:
       AbstractSurfaceLocation
    ts::TS
    u::A
    z::FT
end

"""
  ValuesInterior{FT, A, TS <: TD.ThermodynamicState}

Input container for state variables at the first interior node.

# Fields

$(DSE.FIELDS)
"""
struct ValuesInterior{FT, A, TS <: TD.ThermodynamicState} <:
       AbstractSurfaceLocation
    ts::TS
    u::A
    z::FT
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

qt_in(sc::AbstractSurfaceConditions) = total_specific_humidity(ts_in(sc))
qt_sfc(sc::AbstractSurfaceConditions) = total_specific_humidity(ts_sfc(sc))
Δqt(sc::AbstractSurfaceConditions) = qt_in(sc) - qt_sfc(sc)

b_in(sc::AbstractSurfaceConditions) = sc.state_in.u
b_sfc(sc::AbstractSurfaceConditions) = sc.state_sfc.u

function windspeed(sc::AbstractSurfaceConditions)
    return max(hypot(Δu1(sc), Δu2(sc)), sc.gustiness)
end

"""
  surface_conditions(param_set, sc, scheme)
Function for the computation of surface conditions
based on the Monin-Obukhov similarity functions. Requires
information about thermodynamic parameters (`param_set`)
the surface state `sc`, and the discretisation `scheme`.

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
    param_set::APS,
    sc::AbstractSurfaceConditions,
    scheme::SolverScheme = FVScheme(),
) where {APS, F}
    # Compute output variables for SurfaceFluxConditions
    FT = eltype(sc.state_in.u)
    L_MO = obukhov_length(param_set, sc, scheme)
    ustar = compute_ustar(param_set, L_MO, sc, scheme)
    Cd = momentum_exchange_coefficient(param_set, L_MO, sc, scheme)
    Ch = heat_exchange_coefficient(param_set, L_MO, sc, scheme)
    shf = sensible_heat_flux(param_set, Ch, sc, scheme)
    lhf = latent_heat_flux(param_set, Ch, sc, scheme)
    ρτxz, ρτyz = momentum_fluxes(param_set, Cd, sc, scheme)
    return SurfaceFluxConditions{FT}(L_MO, shf, lhf, ρτxz, ρτyz, ustar, Cd, Ch)
end

function local_lmo(param_set, x_lmo, sc::Union{Fluxes}, scheme)
    κ = CPSGS.von_karman_const(param_set)
    u_scale = compute_physical_scale_coeff(
        param_set,
        sc,
        x_lmo,
        UF.MomentumTransport(),
        scheme,
    )
    return -(windspeed(sc) * u_scale)^3 / κ /
           compute_buoyancy_flux(param_set, sc, scheme)
end

function local_lmo(param_set, x_lmo, sc::ValuesOnly, scheme)
    κ = CPSGS.von_karman_const(param_set)
    u_scale = compute_physical_scale_coeff(
        param_set,
        sc,
        x_lmo,
        UF.MomentumTransport(),
        scheme,
    )
    return (windspeed(sc) * u_scale)^2 / κ /
           compute_bstar(param_set, x_lmo, sc, scheme)
end

"""
  obukhov_length(param_set, sc, scheme)
Compute and return the Monin-Obukhov lengthscale. This is the most 
general case, and assumes non-linear solver iterations are required 
to determine the solution.
"""
function obukhov_length(param_set, sc::AbstractSurfaceConditions, scheme)
    FT = eltype(sc.L_MO_init)
    local sol
    function root_l_mo(x_lmo)
        root = x_lmo - local_lmo(param_set, x_lmo, sc, scheme)
        return root
    end
    sol = RS.find_zero(
        root_l_mo,
        RS.NewtonsMethodAD(sc.L_MO_init),
        RS.CompactSolution(),
        RS.SolutionTolerance(sqrt(eps(FT))),
        100,
    )
    if sol.converged
        L_MO = sol.root
    else
        KA.@print("Warning: Unconverged Surface Fluxes. Returning neutral condition solution \n")
        L_MO = sol.root
    end
    return FT(L_MO)
end

"""
  obukhov_length(param_set, sc, scheme)
Compute and return the Monin-Obukhov lengthscale (LMO) when surface fluxes
and friction velocity are known. This is a case-specific implementation of `obukhov_length`;
iterations to determine LMO are not required.
"""
function obukhov_length(param_set, sc::FluxesAndFrictionVelocity, scheme)
    return -sc.ustar^3 / CPSGS.von_karman_const(param_set) /
           compute_buoyancy_flux(param_set, sc, scheme)
end

"""
  obukhov_length(param_set, sc, scheme)
Compute and return the Monin-Obukhov lengthscale (LMO) when exchange coefficients are known. 
This is a case-specific implementation of `obukhov_length`;
iterations to determine LMO are not required.
"""
function obukhov_length(param_set, sc::Coefficients, scheme)
    FT = eltype(sc.Ch)
    Δθ = TD.dry_pottemp(ts_in(sc)) - TD.dry_pottemp(ts_sfc(sc))
    ρ_sfc = TD.air_density(sc.state_sfc.ts)
    cp_m::FT = TD.cp_m(sc.state_in.ts)
    L_v::FT = TD.latent_heat_vapor(sc.state_in.ts)
    shf = -ρ_sfc * sc.Ch * cp_m * windspeed(sc) * Δθ
    lhf = -ρ_sfc * sc.Ch * L_v * windspeed(sc) * Δqt(sc)
    ustar = sqrt(sc.Cd) * windspeed(sc)
    intermediate_surface_condition = GivenFluxesAndFrictionVelocity(
        sc.state_in,
        sc.state_sfc,
        shf,
        lhf,
        ustar,
        sc.uf,
        sc.z0m,
        sc.z0b,
    )
    return obukhov_length(param_set, intermediate_surface_condition, scheme)
end

function compute_buoyancy_flux(param_set, sc::Coefficients, scheme)
    FT = eltype(sc.Cd)
    # Unpack planet parameters
    grav::FT = CPP.grav(param_set)
    ε_vd::FT = CPP.molmass_ratio(param_set)
    # Compute thermodynamic variables
    cp_m::FT = TD.cp_m(sc.state_in.ts)
    L_v::FT = TD.latent_heat_vapor(sc.state_in.ts)
    # Compute fluxes
    Δθ = TD.dry_pottemp(ts_in(sc)) - TD.dry_pottemp(ts_sfc(sc))
    J = -sc.Ch * cp_m * windspeed(sc) * Δθ
    D = -sc.Ch * L_v * windspeed(sc) * Δqt(sc)
    # Return <w'b'>
    return grav / T_in * (J / cp_m + (ε_vd - 1) * D / L_v)
end
function compute_buoyancy_flux(
    param_set,
    sc::Union{FluxesAndFrictionVelocity, Fluxes, ValuesOnly},
    scheme,
)
    FT = eltype(sc.state_in.u)
    grav::FT = CPP.grav(param_set)
    ε_vd::FT = CPP.molmass_ratio(param_set)

    cp_m::FT = TD.cp_m(sc.state_in.ts)
    L_v::FT = TD.latent_heat_vapor(sc.state_in.ts)
    ρ_sfc = TD.air_density(sc.state_sfc.ts)
    T_in = TD.air_temperature(sc.state_in.ts)

    J = sc.shf / ρ_sfc
    D = sc.lhf / ρ_sfc

    return grav / T_in * (J / cp_m + (ε_vd - 1) * D / L_v)
end

compute_bstar(param_set, L_MO, sc::FluxesAndFrictionVelocity, scheme) =
    -compute_buoyancy_flux(param_set, sc, scheme) / sc.ustar

compute_bstar(param_set, L_MO, sc::Fluxes, scheme) =
    -compute_buoyancy_flux(param_set, sc, scheme) /
    compute_ustar(param_set, L_MO, sc, scheme)

function compute_bstar(param_set, L_MO, sc::AbstractSurfaceConditions, scheme)
    FT = eltype(L_MO)
    uf = sc.uf(param_set, L_MO)
    grav::FT = CPP.grav(param_set)

    ρ_in = TD.air_density(ts_in(sc))
    ρ_sfc = TD.air_density(ts_sfc(sc))

    ρ_star =
        compute_physical_scale_coeff(
            param_set,
            sc,
            L_MO,
            UF.HeatTransport(),
            scheme,
        ) * (ρ_in - ρ_sfc)

    return -grav / ρ_in * ρ_star
end

"""
  compute_ustar(param_set, L_MO, sc, scheme)
Return the friction velocity ustar.
"""
compute_ustar(param_set, L_MO, sc::FluxesAndFrictionVelocity, scheme) = sc.ustar

"""
  compute_ustar(param_set, L_MO, sc, scheme)
Compute and return the friction velocity given the 
Monin-Obukhov lengthscale.
"""
compute_ustar(param_set, L_MO, sc::Fluxes, scheme) =
    windspeed(sc) * compute_physical_scale_coeff(
        param_set,
        sc,
        L_MO,
        UF.MomentumTransport(),
        scheme,
    )

"""
  compute_ustar(param_set, L_MO, sc, scheme)
Compute and return the friction velocity given the 
exchange coefficients.
"""
compute_ustar(param_set, L_MO, sc::Coefficients, scheme) =
    sqrt(sc.Cd) * (windspeed(sc))

"""
  compute_ustar(param_set, L_MO, sc, scheme)
Compute and return the friction velocity given the Monin-Obukhov
lengthscale.
"""
compute_ustar(param_set, L_MO, sc::ValuesOnly, scheme) =
    windspeed(sc) * compute_physical_scale_coeff(
        param_set,
        sc,
        L_MO,
        UF.MomentumTransport(),
        scheme,
    )

"""
  momentum_exchange_coefficient(param_set, L_MO, sc, scheme)
Compute and return Cd, the momentum exchange coefficient, given the 
Monin-Obukhov lengthscale.
"""
function momentum_exchange_coefficient(
    param_set,
    L_MO,
    sc::Union{Fluxes, ValuesOnly, FluxesAndFrictionVelocity},
    scheme,
)
    κ = CPSGS.von_karman_const(param_set)
    transport = UF.MomentumTransport()
    FT = eltype(L_MO)
    ρ_in = TD.air_density(ts_in(sc))
    ρ_sfc = TD.air_density(ts_sfc(sc))
    Δρ = ρ_in - ρ_sfc
    ustar = compute_ustar(param_set, L_MO, sc, scheme)
    if Δρ ≈ FT(0)
        Cd = (κ / log(Δz(sc) / z0(sc, transport)))^2
    else
        Cd = ustar^2 / windspeed(sc)^2
    end
    return Cd
end

"""
  momentum_exchange_coefficient(param_set, L_MO, sc, scheme)
Return Cd, the momentum exchange coefficient.
"""
function momentum_exchange_coefficient(
    param_set,
    L_MO,
    sc::Coefficients,
    scheme,
)
    return sc.Cd
end

"""
  heat_exchange_coefficient(param_set, L_MO, sc, scheme)
Compute and return Ch, the heat exchange coefficient given the 
Monin-Obukhov lengthscale.
"""
function heat_exchange_coefficient(
    param_set,
    L_MO,
    sc::Union{Fluxes, ValuesOnly, FluxesAndFrictionVelocity},
    scheme,
)
    FT = eltype(L_MO)
    Δθ = TD.dry_pottemp(ts_in(sc)) - TD.dry_pottemp(ts_sfc(sc))
    ustar = compute_ustar(param_set, L_MO, sc, scheme)
    θstar =
        Δθ * compute_physical_scale_coeff(
            param_set,
            sc,
            L_MO,
            UF.HeatTransport(),
            scheme,
        )
    if Δθ ≈ FT(0) || windspeed(sc) ≈ FT(0)
        Ch = FT(0)
    else
        Ch = ustar * θstar / windspeed(sc) / Δθ
    end
    return Ch
end

"""
  heat_exchange_coefficient(param_set, L_MO, sc, scheme)
Return Ch, the heat exchange coefficient given the 
Monin-Obukhov lengthscale.
"""
function heat_exchange_coefficient(param_set, L_MO, sc::Coefficients, scheme)
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
  - scheme: Discretisation scheme (currently supports DG and FV)
"""
function momentum_fluxes(param_set, Cd, sc::AbstractSurfaceConditions, scheme)
    ρ_sfc = TD.air_density(ts_sfc(sc))
    ρτxz = -ρ_sfc * Cd * Δu1(sc) * windspeed(sc)
    ρτyz = -ρ_sfc * Cd * Δu2(sc) * windspeed(sc)
    return (ρτxz, ρτyz)
end

"""
  sensible_heat_flux(param_set, Ch, sc, scheme)
In cases where surface fluxes are known, 
return the known latent heat flux.
## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - Ch: Thermal exchange coefficient
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - scheme: Discretisation scheme (currently supports DG and FV)
"""
function sensible_heat_flux(
    param_set,
    Ch,
    sc::Union{ValuesOnly, Coefficients},
    scheme,
)
    FT = eltype(sc.state_in.u)
    grav::FT = CPP.grav(param_set)
    cp_d::FT = CPP.cp_d(param_set)
    R_d::FT = CPP.R_d(param_set)
    T_0::FT = CPP.T_0(param_set)
    cp_m = TD.cp_m(sc.state_in.ts)
    ρ_sfc = TD.air_density(sc.state_sfc.ts)
    T_in = TD.air_temperature(ts_in(sc))
    T_sfc = TD.air_temperature(ts_sfc(sc))
    ΔT = T_in - T_sfc
    hd_sfc = cp_d * (T_sfc - T_0) + R_d * T_0
    ΔΦ = grav * Δz(sc)
    E = -ρ_sfc * Ch * windspeed(sc) * Δqt(sc)
    return -ρ_sfc * Ch * windspeed(sc) * (cp_m * ΔT + ΔΦ) - (hd_sfc) * E
end


"""
  latent_heat_flux(param_set, Ch, sc, scheme)
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
In cases where surface values, or exchange coefficients, are known, compute and
return the latent heat flux. 
"""
function latent_heat_flux(
    param_set,
    Ch,
    sc::Union{ValuesOnly, Coefficients},
    scheme,
)
    FT = eltype(sc.state_in.u)
    grav::FT = CPP.grav(param_set)
    ρ_sfc = TD.air_density(sc.state_in.ts)
    cp_v::FT = CPP.cp_v(param_set)
    Lv_0::FT = CPP.LH_v0(param_set)
    T_in = TD.air_temperature(ts_in(sc))
    T_sfc = TD.air_temperature(ts_sfc(sc))
    ΔT = T_in - T_sfc
    hv_sfc = cp_v * (ΔT) + Lv_0
    Φ_sfc = grav * z_sfc(sc)
    E = -ρ_sfc * Ch * windspeed(sc) * Δqt(sc)
    lhf = (hv_sfc + Φ_sfc) * E
    return lhf
end


"""
  compute_physical_scale_coeff(param_set, sc, L_MO, transport, ::DGScheme)
Computes the coefficient for the physical scale of a variable based on Nishizawa(2018)
for the FV scheme. 

## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination 
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - L_MO: Monin-Obukhov length
  - transport: Transport type, (e.g. Momentum or Heat, used to determine physical scale coefficients)
  - scheme: Discretisation scheme (currently supports DG and FV)
"""
function compute_physical_scale_coeff(
    param_set,
    sc::AbstractSurfaceConditions,
    L_MO::FT,
    transport,
    ::FVScheme,
) where {FT}
    von_karman_const::FT = CPSGS.von_karman_const(param_set)
    uf = sc.uf(param_set, L_MO)
    _π_group = FT(UF.π_group(uf, transport))
    _π_group⁻¹ = (1 / _π_group)
    R_z0 = 1 - z0(sc, transport) / Δz(sc)
    denom1 = log(Δz(sc) / z0(sc, transport))
    denom2 = -UF.Psi(uf, Δz(sc) / uf.L, transport)
    denom3 =
        z0(sc, transport) / Δz(sc) *
        UF.Psi(uf, z0(sc, transport) / uf.L, transport)
    denom4 = R_z0 * (UF.psi(uf, z0(sc, transport) / uf.L, transport) - 1)
    Σterms = denom1 + denom2 + denom3 + denom4
    return _π_group⁻¹ * von_karman_const / Σterms
end


"""
  compute_physical_scale_coeff(param_set, sc, L_MO, transport, ::DGScheme)
Computes the coefficient for the physical scale of a variable based on Byun(1990)
for the DG scheme. 

## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination 
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - L_MO: Monin-Obukhov length
  - transport: Transport type, (e.g. Momentum or Heat, used to determine physical scale coefficients)
  - scheme: Discretisation scheme (currently supports DG and FV)
"""
function compute_physical_scale_coeff(
    param_set,
    sc::AbstractSurfaceConditions,
    L_MO::FT,
    transport,
    ::DGScheme,
) where {FT}
    von_karman_const::FT = CPSGS.von_karman_const(param_set)
    uf = sc.uf(param_set, L_MO)
    _π_group = FT(UF.π_group(uf, transport))
    _π_group⁻¹ = (1 / _π_group)
    denom1 = log(Δz(sc) / z0(sc, transport))
    denom2 = -UF.psi(uf, Δz(sc) / uf.L, transport)
    denom3 = UF.psi(uf, z0(sc, transport) / uf.L, transport)
    Σterms = denom1 + denom2 + denom3
    return _π_group⁻¹ * von_karman_const / Σterms
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
  - X_in,X_sfc: For variable X, values at interior and surface nodes
  - transport: Transport type, (e.g. Momentum or Heat, used to determine physical scale coefficients)
  - scheme: Discretisation scheme (currently supports DG and FV)
"""
function recover_profile(
    param_set::APS,
    sc::AbstractSurfaceConditions,
    L_MO,
    Z,
    X_in,
    X_sfc,
    transport,
    scheme::Union{FVScheme, DGScheme},
) where {FT}
    @assert isless.(Z, sc.vals_in.z)
    uf = sc.uf(param_set, L_MO)
    von_karman_const::FT = CPSGS.von_karman_const(param_set)
    _π_group = FT(UF.π_group(uf, transport))
    _π_group⁻¹ = (1 / _π_group)
    num1 = log(Z / z0(sc, transport))
    num2 = -UF.psi(uf, Z / L_MO, transport)
    num3 = UF.psi(uf, z0(sc, transport) / L_MO, transport)
    Σnum = num1 + num2 + num3
    ΔX = X_in - X_sfc
    return Σnum *
           compute_physical_scale_coeff(
               param_set,
               sc,
               L_MO,
               transport,
               scheme,
           ) *
           _π_group⁻¹ *
           ΔX + X_sfc
end

monin_obukhov_length(sfc::SurfaceFluxConditions) = sfc.L_MO

"""
    get_scalar_flux_coefficients(
       param_set,
       ts_in, ts_sfc,
       z_in, z_sfc,
       u_in, u_sfc,
       L_MO,
       z0m::FT,
       scheme::SolverScheme,
       universal_func,
    )

Returns the exchange coefficients for bulk transfer formulas associated
with the frictional parameters `x_star` and the surface layer similarity
profiles. Taken from equations (36) and (37) of Byun (1990).
"""
function get_scalar_flux_coefficients(
    param_set,
    sc::AbstractSurfaceConditions,
    L_MO,
    z0b::FT,
    scheme::SolverScheme,
) where {VFT, FT}
    N = length(x_in)
    Δx = x_in - x_sfc
    ρ_in = TD.air_density(ts_in)
    ρ_sfc = TD.air_density(ts_sfc)
    Δρ = ρ_in - ρ_sfc # TODO Scalar vs buoyancy / density check for zero-flux condition
    if Δρ ≈ FT(0)
        C = FT(0)
    else
        u_scale = compute_physical_scale_coeff(
            param_set,
            sc,
            L_MO,
            UF.MomentumTransport(),
            scheme,
        )
        x_scale = compute_physical_scale_coeff(
            param_set,
            sc,
            L_MO,
            UF.MomentumTransport(),
            scheme,
        )
        C = u_scale * x_scale
    end
    return C
end

end # SurfaceFluxes module
