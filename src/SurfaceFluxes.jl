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
const TD = Thermodynamics

import StaticArrays
const SA = StaticArrays

include("UniversalFunctions.jl")
import .UniversalFunctions
const UF = UniversalFunctions

include("Parameters.jl")
import .Parameters

using QuadGK

const SFP = Parameters
const APS = SFP.AbstractSurfaceFluxesParameters

abstract type SolverScheme end
struct FVScheme <: SolverScheme end
struct FDScheme <: SolverScheme end


abstract type CanopyType end
struct SparseCanopy{FT} <: CanopyType
    d::FT
    z_star::FT
end
struct DenseCanopy{FT} <: CanopyType
    d::FT
    z_star::FT
end

abstract type AbstractRoughnessSublayerType end
struct NoRSL <: AbstractRoughnessSublayerType end
struct PhysickRSL <: AbstractRoughnessSublayerType
    canopy::CanopyType
end
struct DeRidderRSL <: AbstractRoughnessSublayerType
    canopy::CanopyType
end
struct HarmonRSL <: AbstractRoughnessSublayerType
    canopy::CanopyType
end

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
    SurfaceValues

Input container for state variables at the ground level.

# Fields

$(DSE.FIELDS)
"""
struct SurfaceValues{FT <: Real, A, TS <: TD.ThermodynamicState}
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
struct InteriorValues{FT <: Real, A, TS <: TD.ThermodynamicState}
    z::FT
    u::A
    ts::TS
end

abstract type AbstractSurfaceConditions{FT <: Real, VI <: InteriorValues, VS <: SurfaceValues} end


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
    beta::FT = FT(1)
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
    gustiness::FT = FT(1)
    beta::FT = FT(1)
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

qt_in(param_set::APS, sc::AbstractSurfaceConditions) =
    TD.total_specific_humidity(SFP.thermodynamics_params(param_set), ts_in(sc))
qt_sfc(param_set::APS, sc::AbstractSurfaceConditions) =
    TD.total_specific_humidity(SFP.thermodynamics_params(param_set), ts_sfc(sc))
Δqt(param_set::APS, sc::AbstractSurfaceConditions) = qt_in(param_set, sc) - qt_sfc(param_set, sc)

u_in(sc::AbstractSurfaceConditions) = sc.state_in.u
u_sfc(sc::AbstractSurfaceConditions) = sc.state_sfc.u

function windspeed(sc::AbstractSurfaceConditions)
    return max(hypot(Δu1(sc), Δu2(sc)), sc.gustiness)
end

"""
    surface_conditions(
        param_set::AbstractSurfaceFluxesParameters,
        sc::SurfaceFluxes.AbstractSurfaceConditions{FT},
        scheme::SurfaceFluxes.SolverScheme = FVScheme();
        tol::RS.AbstractTolerance = RS.RelativeOrAbsoluteSolutionTolerance(FT(3e-5), FT(Δz(sc) / 50)),
        tol_neutral::FT = SFP.cp_d(param_set) / 100,
        maxiter::Int = 10,
        soltype::RS.SolutionType = RS.CompactSolution(),
    ) where {FT}

The main user facing function of the module.
It computes the surface conditions
based on the Monin-Obukhov similarity functions. Requires
information about thermodynamic parameters (`param_set`)
the surface state `sc`, the universal function type and
the discretisation `scheme`. Default tolerance for 
Monin-Obukhov length is absolute (i.e. has units [m]).
Returns the RootSolvers `CompactSolution` by default.

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
    sc::AbstractSurfaceConditions{FT},
    scheme::SolverScheme = FVScheme();
    tol::RS.AbstractTolerance = RS.RelativeOrAbsoluteSolutionTolerance(FT(3e-5), FT(Δz(sc) / 50)),
    tol_neutral = FT(SFP.cp_d(param_set) / 100),
    maxiter::Int = 10,
    soltype::RS.SolutionType = RS.CompactSolution(),
    noniterative_stable_sol::Bool = false,
) where {FT}
    uft = SFP.universal_func_type(param_set)
    L_MO = obukhov_length(param_set, sc, uft, scheme; tol, tol_neutral, maxiter, soltype, noniterative_stable_sol)
    ustar = compute_ustar(param_set, L_MO, sc, uft, scheme)
    Cd = momentum_exchange_coefficient(param_set, L_MO, sc, uft, scheme, tol_neutral)
    Ch = heat_exchange_coefficient(param_set, L_MO, sc, uft, scheme, tol_neutral)
    shf = sensible_heat_flux(param_set, Ch, sc, scheme)
    lhf = latent_heat_flux(param_set, Ch, sc, scheme)
    buoy_flux = compute_buoyancy_flux(param_set, shf, lhf, ts_in(sc), ts_sfc(sc), scheme)
    ρτxz, ρτyz = momentum_fluxes(param_set, Cd, sc, scheme)
    E = evaporation(param_set, sc, Ch)
    return SurfaceFluxConditions{FT}(L_MO, shf, lhf, buoy_flux, ρτxz, ρτyz, ustar, Cd, Ch, E)
end

"""
    obukhov_length(sfc::SurfaceFluxConditions)

    obukhov_length( # internal method
        param_set::AbstractSurfaceFluxesParameters,
        sc::AbstractSurfaceConditions,
        uft,
        scheme;
        tol::RS.AbstractTolerance = RS.RelativeOrAbsoluteSolutionTolerance(FT(3e-5), FT(Δz(sc) / 50)),
        tol_neutral::FT = SFP.cp_d(param_set) / 100,
        maxiter::Int = 10
        soltype::RS.SolutionType = RS.CompactSolution(),
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

function non_zero(value::FT) where {FT}
    if abs(value) < eps(FT)
        return value + sqrt(eps(FT))
    else
        return value
    end
end

function obukhov_length(
    param_set,
    sc::AbstractSurfaceConditions{FT},
    uft::UF.AUFT,
    scheme;
    tol::RS.AbstractTolerance = RS.RelativeOrAbsoluteSolutionTolerance(FT(3e-5), FT(Δz(sc) / 50)),
    tol_neutral = FT(SFP.cp_d(param_set) / 100),
    maxiter::Int = 10,
    soltype::RS.SolutionType = RS.CompactSolution(),
    noniterative_stable_sol::Bool = false,
) where {FT}
    thermo_params = SFP.thermodynamics_params(param_set)
    ufparams = SFP.uf_params(param_set)
    grav = SFP.grav(param_set)
    cp_d = SFP.cp_d(param_set)
    DSEᵥ_in = TD.virtual_dry_static_energy(thermo_params, ts_in(sc), grav * z_in(sc))
    DSEᵥ_sfc = TD.virtual_dry_static_energy(thermo_params, ts_sfc(sc), grav * z_sfc(sc))
    ΔDSEᵥ = DSEᵥ_in - DSEᵥ_sfc
    if ΔDSEᵥ >= FT(0) && noniterative_stable_sol == true # Stable Layer
        ### Analytical Solution 
        ### Gryanik et al. (2021)
        ### DOI: 10.1029/2021MS002590)
        Ri_b = (grav * z_in(sc) * ΔDSEᵥ) / (DSEᵥ_sfc * (windspeed(sc))^2)
        @assert Ri_b >= FT(0)
        ϵₘ = FT(Δz(sc) / z0(sc, UF.MomentumTransport()))
        ϵₕ = FT(Δz(sc) / z0(sc, UF.HeatTransport()))
        C = (log(ϵₘ))^2 / (log(ϵₕ))
        ζₐ = ufparams.ζ_a
        ufₐ = UF.universal_func(uft, Δz(sc) / ζₐ, SFP.uf_params(param_set))
        γ = ufparams.γ
        Pr₀ = UF.Pr_0(ufₐ)
        ψ_ma = UF.psi(ufₐ, ζₐ, UF.MomentumTransport())
        ψ_ha = UF.psi(ufₐ, ζₐ, UF.HeatTransport())
        A =
            (log(ϵₘ) - ψ_ma)^(2 * (γ - 1)) / (ζₐ^(γ - 1) * (log(ϵₕ) - ψ_ha)^(γ - 1)) *
            ((log(ϵₘ) - ψ_ma)^2 / (log(ϵₕ) - ψ_ha) - C)
        ζₛ = FT(C) * FT(Ri_b) + FT(A) * FT(Ri_b^γ)
        ufₛ = UF.universal_func(uft, Δz(sc) / ζₛ, SFP.uf_params(param_set))
        ψₘ = UF.psi(ufₛ, ζₛ, UF.MomentumTransport())
        ψₕ = UF.psi(ufₛ, ζₛ, UF.HeatTransport())
        # Compute exchange coefficients
        κ = SFP.von_karman_const(param_set)
        Cd = κ^2 / (log(ϵₘ)^2) * (1 - ψₘ / log(ϵₘ))^(-2)
        Ch = κ^2 / (Pr₀ * log(ϵₘ) * log(ϵₕ)) * (1 - ψₘ / log(ϵₘ))^(-1) * (1 - ψₕ / Pr₀ / log(ϵₕ))^(-1)
        return non_zero(Δz(sc) ./ ζₛ)
    elseif tol_neutral >= abs(ΔDSEᵥ) # Neutral Layer
        # Large L_MO -> virtual dry static energy suggests neutral boundary layer
        # Return ζ->0 in the neutral boundary layer case, where ζ = z / L_MO
        return L_MO = FT(Inf * sign(non_zero(ΔDSEᵥ)))
    else
        function root_l_mo(x_lmo)
            residual = x_lmo - local_lmo(param_set, x_lmo, sc, uft, scheme)
            return residual
        end
        L_MO_init = FT(-1)
        sol = RS.find_zero(root_l_mo, RS.NewtonsMethodAD(L_MO_init), soltype, tol, maxiter)
        L_MO = sol.root
        if !sol.converged
            if error_on_non_convergence()
                KA.@print("maxiter reached in SurfaceFluxes.jl:\n")
                KA.@print(" T_in = ", TD.air_temperature(thermo_params, ts_in(sc)))
                KA.@print(", T_sfc = ", TD.air_temperature(thermo_params, ts_sfc(sc)))
                KA.@print(", q_in = ", TD.total_specific_humidity(thermo_params, ts_in(sc)))
                KA.@print(", q_sfc = ", TD.total_specific_humidity(thermo_params, ts_sfc(sc)))
                KA.@print(", u_in = ", u_in(sc))
                KA.@print(", u_sfc = ", u_sfc(sc))
                KA.@print(", z0_m = ", z0(sc, UF.MomentumTransport()))
                KA.@print(", z0_b = ", z0(sc, UF.HeatTransport()))
                KA.@print(", Δz = ", Δz(sc))
                KA.@print(", ΔDSEᵥ = ", ΔDSEᵥ)
                KA.@print("\n")
                KA.@print("ts_in(sc) = ", ts_in(sc))
                KA.@print("\n")
                KA.@print("ts_sfc(sc) = ", ts_sfc(sc))
                KA.@print("\n")
                if soltype isa RS.CompactSolution
                    KA.@print(", sol.root = ", sol.root)
                else
                    KA.@print(", sol.root_history = ", sol.root_history)
                    KA.@print(", sol.err_history = ", sol.err_history)
                end
                error("Unconverged Surface Fluxes.")
            else
                KA.@print("Warning: Unconverged Surface Fluxes. Returning last interation.")
            end
        end
        return non_zero(L_MO)
    end
end

function obukhov_length(param_set, sc::FluxesAndFrictionVelocity{FT}, uft::UF.AUFT, scheme; kwargs...) where {FT}
    return -sc.ustar^3 / FT(SFP.von_karman_const(param_set)) / non_zero(compute_buoyancy_flux(param_set, sc, scheme))
end

function obukhov_length(param_set, sc::Coefficients{FT}, uft::UF.AUFT, scheme; kwargs...) where {FT}
    lhf = latent_heat_flux(param_set, sc.Ch, sc, scheme)
    shf = sensible_heat_flux(param_set, sc.Ch, sc, scheme)
    ustar = sqrt(sc.Cd) * windspeed(sc)
    buoyancy_flux = compute_buoyancy_flux(param_set, shf, lhf, ts_in(sc), ts_sfc(sc), scheme)
    return -ustar^3 / FT(SFP.von_karman_const(param_set)) / non_zero(buoyancy_flux)
end

"""
    local_lmo(param_set, x_lmo, sc, uft, scheme)

Helper function for the iterative solver with the Monin-Obukhov length equation with buoyancy flux.
"""
function local_lmo(param_set, x_lmo, sc::Fluxes, uft::UF.AUFT, scheme)
    κ = SFP.von_karman_const(param_set)
    u_scale = compute_physical_scale_coeff(param_set, sc, x_lmo, UF.MomentumTransport(), uft, scheme)
    return -(windspeed(sc) * u_scale)^3 / κ / non_zero(compute_buoyancy_flux(param_set, sc, scheme))
end

"""
    local_lmo(param_set, x_lmo, sc, uft, scheme)

Helper function for the iterative solver with the Monin-Obukhov
length equation with buoyancy star.
"""
function local_lmo(param_set, x_lmo, sc::ValuesOnly, uft::UF.AUFT, scheme)
    κ = SFP.von_karman_const(param_set)
    u_scale = compute_physical_scale_coeff(param_set, sc, x_lmo, UF.MomentumTransport(), uft, scheme)
    return (windspeed(sc) * u_scale)^2 / κ / non_zero(compute_bstar(param_set, x_lmo, sc, uft, scheme))
end

"""
    compute_buoyancy_flux(param_set, shf, lhf, ts_in, ts_sfc, scheme)

Returns the buoyancy flux when the surface fluxes are known.
"""
function compute_buoyancy_flux(param_set, shf::FT, lhf::FT, ts_in, ts_sfc, scheme) where {FT}
    thermo_params = SFP.thermodynamics_params(param_set)
    grav = SFP.grav(param_set)
    ε_vd = SFP.molmass_ratio(param_set)
    cp_m = TD.cp_m(thermo_params, ts_in)
    L_v = TD.latent_heat_vapor(thermo_params, ts_in)
    ρ_sfc = TD.air_density(thermo_params, ts_sfc)
    T_in = TD.air_temperature(thermo_params, ts_in)
    return grav / ρ_sfc * (shf / cp_m / T_in + (ε_vd - 1) * lhf / L_v)
end

function compute_buoyancy_flux(param_set, sc::Union{FluxesAndFrictionVelocity, Fluxes}, scheme)
    return compute_buoyancy_flux(param_set, sc.shf, sc.lhf, ts_in(sc), ts_sfc(sc), scheme)
end

"""
    compute_bstar(param_set, L_MO, sc, uft, scheme)

Returns buoyancy star based on known air densities.
"""
function compute_bstar(param_set, L_MO, sc::AbstractSurfaceConditions{FT}, uft::UF.AUFT, scheme) where {FT}
    thermo_params = SFP.thermodynamics_params(param_set)
    uf = UF.universal_func(uft, L_MO, SFP.uf_params(param_set))
    grav::FT = SFP.grav(param_set)
    DSEᵥ_sfc = TD.virtual_dry_static_energy(thermo_params, ts_sfc(sc), grav * z_sfc(sc))
    DSEᵥ_in = TD.virtual_dry_static_energy(thermo_params, ts_in(sc), grav * z_in(sc))
    DSEᵥ_star =
        compute_physical_scale_coeff(param_set, sc, L_MO, UF.HeatTransport(), uft, scheme) * (DSEᵥ_in - DSEᵥ_sfc)
    return grav * DSEᵥ_star / DSEᵥ_in
end

"""
    compute_bstar(param_set, L_MO, sc, uft, scheme)

Returns buoyancy star based on known friction velocity  and fluxes.
"""
compute_bstar(param_set, L_MO, sc::FluxesAndFrictionVelocity, uft, scheme) =
    -compute_buoyancy_flux(param_set, sc, scheme) / sc.ustar

"""
    compute_bstar(param_set, L_MO, sc, uft, scheme)

Return buoyancy star based on known fluxes.
"""
compute_bstar(param_set, L_MO, sc::Fluxes, uft, scheme) =
    -compute_buoyancy_flux(param_set, sc, scheme) / compute_ustar(param_set, L_MO, sc, uft, scheme)


"""
    compute_ustar(
        param_set::AbstractSurfaceFluxesParameters,
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
"""
function compute_ustar end

compute_ustar(param_set, L_MO, sc::FluxesAndFrictionVelocity, uft, scheme) = sc.ustar

compute_ustar(param_set, L_MO, sc::Fluxes, uft, scheme) =
    windspeed(sc) * compute_physical_scale_coeff(param_set, sc, L_MO, UF.MomentumTransport(), uft, scheme)

compute_ustar(param_set, L_MO, sc::Coefficients, uft, scheme) = sqrt(sc.Cd) * (windspeed(sc))

function compute_ustar(param_set, L_MO, sc::ValuesOnly, uft, scheme)
    return windspeed(sc) * compute_physical_scale_coeff(param_set, sc, L_MO, UF.MomentumTransport(), uft, scheme)
end

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
    tol_neutral,
)
    FT = eltype(L_MO)
    thermo_params = SFP.thermodynamics_params(param_set)
    κ = SFP.von_karman_const(param_set)
    grav::FT = SFP.grav(param_set)
    transport = UF.MomentumTransport()
    DSEᵥ_sfc = TD.virtual_dry_static_energy(thermo_params, ts_sfc(sc), grav * z_sfc(sc))
    DSEᵥ_in = TD.virtual_dry_static_energy(thermo_params, ts_in(sc), grav * z_in(sc))
    ΔDSEᵥ = DSEᵥ_in - DSEᵥ_sfc
    cp_d::FT = SFP.cp_d(param_set)
    if isapprox(ΔDSEᵥ, FT(0); atol = tol_neutral)
        Cd = (κ / log(Δz(sc) / z0(sc, transport)))^2
    else
        ustar = compute_ustar(param_set, L_MO, sc, uft, scheme)
        Cd = ustar^2 / windspeed(sc)^2
    end
    return Cd
end

"""
    momentum_exchange_coefficient(param_set, L_MO, sc, uft, scheme, tol_neutral)

Return Cd, the momentum exchange coefficient.
"""
function momentum_exchange_coefficient(param_set, L_MO, sc::Coefficients, uft, scheme, tol_neutral)
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
    sc::Union{Fluxes, ValuesOnly, FluxesAndFrictionVelocity},
    uft,
    scheme,
    tol_neutral,
)
    FT = eltype(L_MO)
    thermo_params = SFP.thermodynamics_params(param_set)
    transport = UF.HeatTransport()
    κ = SFP.von_karman_const(param_set)
    grav::FT = SFP.grav(param_set)
    cp_d::FT = SFP.cp_d(param_set)
    DSEᵥ_sfc = TD.virtual_dry_static_energy(thermo_params, ts_sfc(sc), grav * z_sfc(sc))
    DSEᵥ_in = TD.virtual_dry_static_energy(thermo_params, ts_in(sc), grav * z_in(sc))
    ΔDSEᵥ = DSEᵥ_in - DSEᵥ_sfc
    z0_b = z0(sc, UF.HeatTransport())
    z0_m = z0(sc, UF.MomentumTransport())
    if isapprox(ΔDSEᵥ, FT(0); atol = tol_neutral)
        Ch = κ^2 / (log(Δz(sc) / z0_b) * log(Δz(sc) / z0_m))
    else
        ϕ_heat = compute_physical_scale_coeff(param_set, sc, L_MO, transport, uft, scheme)
        ustar = compute_ustar(param_set, L_MO, sc, uft, scheme)
        Ch = ustar * ϕ_heat / windspeed(sc)
    end
    return Ch
end

"""
    heat_exchange_coefficient(param_set, L_MO, sc, uft, scheme)

Return Ch, the heat exchange coefficient given the
Monin-Obukhov lengthscale.
"""
function heat_exchange_coefficient(param_set, L_MO, sc::Coefficients, uft, scheme, tol_neutral)
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
function sensible_heat_flux(param_set, Ch::FT, sc::Union{ValuesOnly, Coefficients}, scheme) where {FT}
    thermo_params = SFP.thermodynamics_params(param_set)
    grav::FT = SFP.grav(param_set)
    cp_d::FT = SFP.cp_d(param_set)
    R_d::FT = SFP.R_d(param_set)
    T_0::FT = SFP.T_0(param_set)
    cp_m = TD.cp_m(thermo_params, ts_in(sc))
    ρ_sfc = TD.air_density(thermo_params, ts_sfc(sc))
    T_in = TD.air_temperature(thermo_params, ts_in(sc))
    T_sfc = TD.air_temperature(thermo_params, ts_sfc(sc))
    ΔT = T_in - T_sfc
    hd_sfc = cp_d * (T_sfc - T_0) + R_d * T_0
    ΔΦ = grav * Δz(sc)
    E = evaporation(param_set, sc, Ch)
    return -ρ_sfc * Ch * windspeed(sc) * (cp_m * ΔT + ΔΦ) - (hd_sfc) * E
end

"""
    sensible_heat_flux(param_set, Ch, sc, scheme)

In cases where surface fluxes are known,
return the known sensible heat flux.
"""
function sensible_heat_flux(param_set, Ch, sc::Union{Fluxes, FluxesAndFrictionVelocity}, scheme)
    return sc.shf
end

"""
    latent_heat_flux(param_set, Ch, sc, scheme)

In cases where surface fluxes are known,
return the known latent heat flux.
"""
function latent_heat_flux(param_set, L_MO, sc::Union{Fluxes, FluxesAndFrictionVelocity}, scheme)
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
function latent_heat_flux(param_set, Ch::FT, sc::Union{ValuesOnly, Coefficients}, scheme) where {FT}
    thermo_params = SFP.thermodynamics_params(param_set)
    grav::FT = SFP.grav(param_set)
    ρ_sfc = TD.air_density(thermo_params, ts_sfc(sc))
    cp_v::FT = SFP.cp_v(param_set)
    Lv_0::FT = SFP.LH_v0(param_set)
    T_0::FT = SFP.T_0(param_set)
    T_sfc = TD.air_temperature(thermo_params, ts_sfc(sc))
    hv_sfc = cp_v * (T_sfc - T_0) + Lv_0
    Φ_sfc = grav * z_sfc(sc)
    E = evaporation(param_set, sc, Ch)
    lhf = (hv_sfc + Φ_sfc) * E
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
function evaporation(param_set, sc::Union{Fluxes, FluxesAndFrictionVelocity}, Ch::FT) where {FT}
    thermo_params = SFP.thermodynamics_params(param_set)
    grav::FT = SFP.grav(param_set)
    T_sfc = TD.air_temperature(thermo_params, ts_sfc(sc))
    hv_sfc = TD.latent_heat_vapor(thermo_params, T_sfc)
    Φ_sfc = grav * z_sfc(sc)
    return sc.lhf / (hv_sfc + Φ_sfc)
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
    param_set::APS,
    sc::AbstractSurfaceConditions,
    L_MO::FT,
    transport,
    uft,
    ::FVScheme,
) where {FT}
    von_karman_const::FT = SFP.von_karman_const(param_set)
    uf = UF.universal_func(uft, L_MO, SFP.uf_params(param_set))
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
  - uft: A Universal Function type, (returned by, e.g. Businger())
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
    von_karman_const::FT = SFP.von_karman_const(param_set)
    uf = UF.universal_func(uft, L_MO, SFP.uf_params(param_set))
    _π_group = FT(UF.π_group(uf, transport))
    _π_group⁻¹ = (1 / _π_group)
    denom1 = log(Δz(sc) / z0(sc, transport))
    denom2 = -UF.psi(uf, Δz(sc) / uf.L, transport)
    denom3 = UF.psi(uf, z0(sc, transport) / uf.L, transport)
    Σterms = denom1 + denom2 + denom3
    return _π_group⁻¹ * von_karman_const / Σterms
end

"""
    recover_profile(param_set, sc, L_MO, Z, X_sfc, transport, uft, scheme, rsl)

Recover profiles of variable X given values of Z coordinates. Follows Nishizawa equation (21,22)
## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - L_MO: Monin-Obukhov length
  - Z: Z coordinate(s) (within surface layer) for which variable values are required
  - X_sfc: For variable X, values at surface nodes
  - transport: Transport type, (e.g. Momentum or Heat, used to determine physical scale coefficients)
  - uft: A Universal Function type, (returned by, e.g., Businger())
  - scheme: Discretization scheme (currently supports FD and FV)
  - rsl : Roughness Sublayer Formulation (e.g. NoRSL, PhysickRSL, DeRidderRSL)
  # TODO: add tests
  # TODO: Verify that all current RSL models fall into this general code pattern, then 𝜙 can be abstracted
"""
function recover_profile(
    param_set::APS,
    sc::AbstractSurfaceConditions,
    L_MO::FT,
    Z,
    X_sfc,
    X_star,
    transport,
    uft::UF.AUFT,
    scheme::Union{FVScheme, FDScheme},
    rsl::NoRSL,
) where {FT}
    uf = UF.universal_func(uft, L_MO, SFP.uf_params(param_set))
    von_karman_const::FT = SFP.von_karman_const(param_set)
    _π_group = FT(UF.π_group(uf, transport))
    _π_group⁻¹ = (1 / _π_group)
    num1 = log(Z / z0(sc, transport))
    num2 = -UF.psi(uf, Z / L_MO, transport)
    num3 = UF.psi(uf, z0(sc, transport) / L_MO, transport)
    Σnum = num1 + num2 + num3

    return Σnum * X_star / von_karman_const + X_sfc
end


"""
    recover_profile(param_set, sc, L_MO, Z, X_sfc, transport, uft, scheme, rsl) 

Recover profiles of variable X given values of Z coordinates, as long as Z > d. 
Follows Nishizawa equation (21, 22) 
Canopy modification follows Physick and Garratt (1995) equation (8, 9)
## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - L_MO: Monin-Obukhov length
  - Z: Z coordinate(s) (within surface layer) for which variable values are required
  - X_sfc: For variable X, values at surface nodes
  - transport: Transport type, (e.g. Momentum or Heat, used to determine physical scale coefficients)
  - uft: A Universal Function type, (returned by, e.g., Businger())
  - scheme: Discretization scheme (currently supports FD and FV)
  - rsl : Roughness Sublayer Formulation (e.g. NoRSL, PhysickRSL, DeRidderRSL)

# TODO: add tests
# TODO: Verify that all current RSL models fall into this general code pattern, then 𝜙 can be abstracted
"""
function recover_profile(
    param_set::APS,
    sc::AbstractSurfaceConditions,
    L_MO::FT,
    Z,
    X_sfc,
    X_star,
    transport,
    uft::UF.AUFT,
    scheme::Union{FVScheme, FDScheme},
    rsl::PhysickRSL,
) where {FT}
    z_star = rsl.canopy.z_star
    d = rsl.canopy.d
    uf = UF.universal_func(uft, L_MO, SFP.uf_params(param_set))
    von_karman_const::FT = SFP.von_karman_const(param_set)
    _π_group = FT(UF.π_group(uf, transport))
    _π_group⁻¹ = (1 / _π_group)
    num1 = log((Z - d) / z0(sc, transport))
    num2 = -UF.psi(uf, (Z - d) / L_MO, transport)
    num3 = UF.psi(uf, z0(sc, transport) / L_MO, transport)
    Σnum = num1 + num2 + num3
    ### Protoype : QuadGK integration to evaluate the canopy RSL effect given by Physick and Garratt (1995)
    ### PG95 assume the same 𝜙 function for momentum and scalar (velocity, heat, moisture etc.) variables
    ### For a model level `Z`, we have the offset coordinate `z = Z-d` over which the functions in PG95 are defined

    # Physick 1995 RSL Term 
    function ψ_RSL_term(z)
        if (z < d || z > z_star)
            return 0
        end

        function integrand(x)
            ϕ = UF.phi(uf, (x - d) / L_MO, transport)

            # Physick1995 Eq. 18
            ϕ_italic = FT(0.5) * exp(log(2) * (x - d) / (z_star - d))
            return ϕ * (1 - ϕ_italic) / (x - d)
        end

        integral, error = quadgk(integrand, z, z_star)
        return integral
    end
    rsl_physick = ψ_RSL_term(Z)
    Σnum += rsl_physick

    return Σnum * X_star / von_karman_const + X_sfc
end

"""
    recover_profile(param_set, sc, L_MO, Z, X_sfc, transport, uft, scheme, rsl) 

Recover profiles of variable X given values of Z coordinates, as long as Z > d. 
Follows Nishizawa equation (21, 22) 
Canopy modification follows De Ridder (2010) equation (6, 7, 13)

## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - L_MO: Monin-Obukhov length
  - Z: Z coordinate(s) (within surface layer) for which variable values are required
  - X_sfc: For variable X, values at surface nodes
  - transport: Transport type, (e.g. Momentum or Heat, used to determine physical scale coefficients)
  - uft: A Universal Function type, (returned by, e.g., Businger())
  - scheme: Discretization scheme (currently supports FD and FV)
  - rsl : Roughness Sublayer Formulation (e.g. NoRSL, PhysickRSL, DeRidderRSL)

# TODO: add tests
# TODO: Verify that all current RSL models fall into this general code pattern, then 𝜙 can be abstracted
"""
function recover_profile(
    param_set::APS,
    sc::AbstractSurfaceConditions,
    L_MO::FT,
    Z,
    X_sfc,
    X_star,
    transport,
    uft::UF.AUFT,
    scheme::Union{FVScheme, FDScheme},
    rsl::DeRidderRSL,
) where {FT}
    z_star = rsl.canopy.z_star
    d = rsl.canopy.d
    uf = UF.universal_func(uft, L_MO, SFP.uf_params(param_set))
    von_karman_const::FT = SFP.von_karman_const(param_set)
    _π_group = FT(UF.π_group(uf, transport))
    _π_group⁻¹ = (1 / _π_group)
    num1 = log((Z - d) / z0(sc, transport))
    num2 = -UF.psi(uf, (Z - d) / L_MO, transport)
    num3 = UF.psi(uf, z0(sc, transport) / L_MO, transport)
    Σnum = num1 + num2 + num3

    ## De Ridder (2010) RSL Term 
    function ψ_analytic(z)
        # TODO Move to ClimaParameters
        # De Ridder Eq. 13

        ν = FT(0.5)
        if transport isa UF.MomentumTransport
            μ = FT(2.59) # Momentum, μₕ = 0.95
        elseif transport isa UF.HeatTransport
            μ = FT(0.95) # Momentum, μₕ = 0.95
        end
        λ = FT(1.5)
        return UF.phi(uf, z / L_MO * (1 + ν / (μ * z / z_star)), transport) *
               (1 / λ) *
               log(1 + (λ / (μ * z / z_star))) *
               exp(-μ * z / z_star)
    end
    rsl_deridder = ψ_analytic(Z - d)
    Σnum += rsl_deridder

    return Σnum * X_star / von_karman_const + X_sfc
end

end # SurfaceFluxes module
