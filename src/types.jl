
"""
    Surface flux configuration specs
"""

abstract type AbstractRoughnessParams end
abstract type AbstractGustinessSpec end



struct ConstantGustinessSpec{TG <: Real} <: AbstractGustinessSpec
    value::TG
end

struct DeardorffGustinessSpec <: AbstractGustinessSpec end

Base.broadcastable(p::AbstractRoughnessParams) = tuple(p)
Base.broadcastable(p::AbstractGustinessSpec) = tuple(p)


struct SurfaceFluxConfig{R <: AbstractRoughnessParams, G <: AbstractGustinessSpec}
    roughness::R
    gustiness::G
end



const FluxOption{FT} = Union{Nothing, FT}

@inline maybe_convert_option(::Nothing, ::Type{FT}) where {FT} = nothing
@inline function maybe_convert_option(value::Number, ::Type{FT}) where {FT}
    return convert(FT, value)
end
@inline function maybe_convert_option(value::FT, ::Type{FT}) where {FT}
    return value
end

struct FluxSpecs{FT}
    shf::FluxOption{FT}
    lhf::FluxOption{FT}
    ustar::FluxOption{FT}
    Cd::FluxOption{FT}
    Ch::FluxOption{FT}
end

FluxSpecs{FT}() where {FT} = FluxSpecs{FT}(nothing, nothing, nothing, nothing, nothing)

function FluxSpecs(::Type{FT};
    shf = nothing,
    lhf = nothing,
    ustar = nothing,
    Cd = nothing,
    Ch = nothing,
) where {FT}
    return FluxSpecs{FT}(
        maybe_convert_option(shf, FT),
        maybe_convert_option(lhf, FT),
        maybe_convert_option(ustar, FT),
        maybe_convert_option(Cd, FT),
        maybe_convert_option(Ch, FT),
    )
end

struct SolverOptions{FT}
    tol::FT
    maxiter::Int
end

function SolverOptions(::Type{FT};
    tol = sqrt(eps(FT)),
    maxiter::Int = 100,
) where {FT}
    return SolverOptions{FT}(convert(FT, tol), maxiter)
end

"""
    SurfaceFluxInputs

Immutable container describing the atmospheric and surface state using primitive
quantities plus module-defined parameterizations. Instances of this type are
passed to the functional surface flux solver.

- `Tin`, `qin`, `ρin`: Interior air temperature [K], specific humidity [kg/kg], and density [kg/m³]
- `Ts_guess`, `qs_guess`: Scalar initial guesses for surface temperature and humidity
- `Φs`: Surface geopotential [m²/s²]
- `Δz`: Height difference between interior and surface reference levels [m]
- `d`: Displacement height [m]
- `u_int`, `u_sfc`: Horizontal wind components (u, v) at interior and surface levels [m/s]
- `roughness_model`: Module-defined roughness parameterization
- `gustiness_model`: Module-defined gustiness parameterization
- `update_Ts!`, `update_qs!`: Optional hooks invoked each solver iteration
- `shf`, `lhf`, `ustar`, `Cd`, `Ch`: Optional prescribed flux/scale quantities supplied via `FluxSpecs`
"""
struct SurfaceFluxInputs{
    FT,
    RM <: AbstractRoughnessParams,
    GM <: AbstractGustinessSpec,
    RI,
    UpdateTs,
    UpdateQs,
    U,
}
    Tin::FT
    qin::FT
    ql_in::FT
    qi_in::FT
    ρin::FT
    Ts_guess::FT
    qs_guess::FT
    Φs::FT
    Δz::FT
    d::FT
    u_int::U
    u_sfc::U
    roughness_model::RM
    gustiness_model::GM
    roughness_inputs::RI
    update_Ts!::UpdateTs
    update_qs!::UpdateQs
    shf::Union{Nothing, FT}
    lhf::Union{Nothing, FT}
    ustar::Union{Nothing, FT}
    Cd::Union{Nothing, FT}
    Ch::Union{Nothing, FT}
end

Base.@kwdef mutable struct SurfaceFluxIterationState{FT}
    Ts::FT = FT(0)
    qs::FT = FTSurfaceFluxIterationState(0)
    ustar::FT = FT(0.1)
    Cd::FT = FT(0)
    Ch::FT = FT(0)
    ρ_sfc::FT = FT(1)
    buoyancy_flux::FT = FT(0)
end

struct CallableContext{FT, U}
    Tin::FT
    qin::FT
    ρin::FT
    Ts::FT
    qs::FT
    Φs::FT
    Δz::FT
    d::FT
    u_int::U
    u_sfc::U
    buoyancy_flux::FT
    ustar::FT
    Cd::FT
    Ch::FT
    ρ_sfc::FT
end

Base.propertynames(::CallableContext) = (
    :Tin,
    :qin,
    :ρin,
    :Ts,
    :qs,
    :Φs,
    :Δz,
    :d,
    :u_int,
    :u_sfc,
    :buoyancy_flux,
    :ustar,
    :Cd,
    :Ch,
    :ρ_sfc,
)

"""
    SurfaceFluxConditions

Surface flux conditions, returned from `surface_fluxes`.

# Fields
- `L_MO::FT`: Monin-Obukhov lengthscale [m]
- `shf::FT`: Sensible heat flux [W/m²]
- `lhf::FT`: Latent heat flux [W/m²]
- `buoy_flux::FT`: Buoyancy flux [m²/s³]
- `ρτxz::FT`: Momentum flux, eastward component [kg/(m·s²)]
- `ρτyz::FT`: Momentum flux, northward component [kg/(m·s²)]
- `ustar::FT`: Friction velocity [m/s]
- `Cd::FT`: Momentum exchange coefficient
- `Ch::FT`: Heat exchange coefficient
- `evaporation::FT`: Evaporation rate [kg/(m²·s)]
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

struct StateValues{FT, TS}
    z::FT
    u::Tuple{FT, FT}
    ts::TS
end

struct Fluxes{S, FT}
    state_int::S
    state_sfc::S
    shf::FT
    lhf::FT
    z0m::FT
    z0h::FT
end

struct FluxesAndFrictionVelocity{S, FT}
    state_int::S
    state_sfc::S
    shf::FT
    lhf::FT
    ustar::FT
    z0m::FT
    z0h::FT
end

struct ValuesOnly{S, FT}
    state_int::S
    state_sfc::S
    z0m::FT
    z0h::FT
end

struct Coefficients{S, FT}
    state_int::S
    state_sfc::S
    Cd::FT
    Ch::FT
    z0m::FT
    z0h::FT
end
