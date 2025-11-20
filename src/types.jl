
### Stability Correction Function Approximations
# 1) `LayerAverageScheme`: Follows Nishizawa & Kitamura (2018) for FV approximations
# 2) `PointValueScheme`: Standard finite difference stencil assumption
abstract type SolverScheme end
struct LayerAverageScheme <: SolverScheme end
struct PointValueScheme <: SolverScheme end

### Roughness Models
# 1) ScalarRoughness (User prescribed constant values)
# 2) CharnockRoughness (Charnock u★ dependent formulation)
abstract type RoughnessModel end
struct CharnockRoughness <: RoughnessModel end
struct ScalarRoughness <: RoughnessModel end
struct FunctionalRoughness <: RoughnessModel end

### Surface Temperature Models
# 1) ScalarTemperature (User prescribed values)
# 2) CharnockRoughness (Charnock u★ dependent formulation)
abstract type SurfaceTemperatureModel end
struct ScalarTemperature <: SurfaceTemperatureModel end
struct DynamicSurfaceTemperature <: SurfaceTemperatureModel end


### Input Variable Containers

"""
   StateValues

Input container for state variables at either first / interior nodes.

# Fields

$(DSE.FIELDS)
"""
struct StateValues{FT <: Real, A, TS <: TD.ThermodynamicState, NT}
    z::FT
    u::A
    ts::TS
    args::NT
end
function StateValues(z::FT, u::A, ts::TS; args::NT = nothing) where {FT, A, TS, NT}
    return StateValues{FT, A, TS, NT}(z,u,ts,args)
end

### Input Containers for surface condtions
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
        roughness_model,
    )
end


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
