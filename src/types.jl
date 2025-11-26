
### Stability Correction Function Approximations
# 1) `LayerAverageScheme`: Follows Nishizawa & Kitamura (2018) for FV approximations
# 2) `PointValueScheme`: Standard finite difference stencil assumption
abstract type SolverScheme end
struct LayerAverageScheme <: SolverScheme end
struct PointValueScheme <: SolverScheme end

"""
    Gustiness and quantity specifications
"""
abstract type GustinessModel end

struct ConstantGustiness{FT} <: GustinessModel
    value::FT
end

struct FunctionalGustiness{F} <: GustinessModel
    fn::F
end

abstract type SurfaceQuantity end

struct SurfaceScalar{FT} <: SurfaceQuantity
    value::FT
end

struct SurfaceCallable{F} <: SurfaceQuantity
    fn::F
end

struct SurfaceSpec{Kind, T}
    value::T
end

surface_temperature(spec) = SurfaceSpec{:temperature, typeof(spec)}(spec)
surface_specific_humidity(spec) = SurfaceSpec{:specific_humidity, typeof(spec)}(spec)
momentum_roughness(spec) = SurfaceSpec{:momentum_roughness, typeof(spec)}(spec)
scalar_roughness(spec) = SurfaceSpec{:scalar_roughness, typeof(spec)}(spec)

@inline surface_quantity(spec::SurfaceSpec{Kind, T}, ::Type{FT}) where {Kind, T, FT} =
    surface_quantity(spec.value, FT)

struct CharnockMomentum{FT}
    α::FT
    grav::FT
end

@inline function (cm::CharnockMomentum{FT})(ctx) where {FT}
    return cm.α * ctx.ustar^2 / cm.grav
end

@inline surface_quantity(x::Number, ::Type{FT}) where {FT} =
    SurfaceScalar(convert(FT, x))
@inline function surface_quantity(spec::SurfaceScalar, ::Type{FT}) where {FT}
    return SurfaceScalar(convert(FT, spec.value))
end
@inline function surface_quantity(spec::SurfaceCallable, ::Type{FT}) where {FT}
    return spec
end
@inline surface_quantity(fn::F, ::Type{FT}) where {F, FT} = SurfaceCallable(fn)

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
    tol_neutral::FT
    maxiter::Int
end

function SolverOptions(::Type{FT};
    tol = sqrt(eps(FT)),
    tol_neutral = sqrt(eps(FT)),
    maxiter::Int = 30,
) where {FT}
    return SolverOptions{FT}(convert(FT, tol), convert(FT, tol_neutral), maxiter)
end

"""
    SurfaceFluxInputs

Immutable container describing the atmospheric and surface state using primitive
quantities or module-defined callables. Instances of this type are passed to the
functional surface flux solver.

- `Tin`, `qin`, `ρin`: Interior air temperature [K], specific humidity [kg/kg], and density [kg/m³]
- `Ts_spec`, `qs_spec`: Surface temperature and specific humidity specifications (scalar or callable)
- `Φs`: Surface geopotential [m²/s²]
- `Δz`: Height difference between interior and surface reference levels [m]
- `d`: Displacement height [m]
- `u_in`, `u_sfc`: Horizontal wind components (u, v) at interior and surface levels [m/s]
- `gustiness_model`: Module-defined gustiness model (constant or callable)
- `z0m_spec`, `z0b_spec`: Momentum and scalar roughness specifications
- `shf`, `lhf`, `ustar`, `Cd`, `Ch`: Optional prescribed flux/scale quantities
- Optional flux/scale quantities are supplied via `FluxSpecs`
"""
struct SurfaceFluxInputs{
    FT,
    TsSpec <: SurfaceQuantity,
    QsSpec <: SurfaceQuantity,
    Z0mSpec <: SurfaceQuantity,
    Z0bSpec <: SurfaceQuantity,
    GM <: GustinessModel,
    U,
}
    Tin::FT
    qin::FT
    ρin::FT
    Ts_spec::TsSpec
    qs_spec::QsSpec
    Φs::FT
    Δz::FT
    d::FT
    u_in::U
    u_sfc::U
    gustiness_model::GM
    z0m_spec::Z0mSpec
    z0b_spec::Z0bSpec
    shf::Union{Nothing, FT}
    lhf::Union{Nothing, FT}
    ustar::Union{Nothing, FT}
    Cd::Union{Nothing, FT}
    Ch::Union{Nothing, FT}
end

function SurfaceFluxInputs(
    Tin::FT,
    qin::FT,
    ρin::FT,
    Ts,
    qs,
    Φs::FT,
    Δz::FT,
    d::FT,
    u_in,
    u_sfc,
    gustiness,
    z0m,
    z0b,
    flux_specs::FluxSpecs{FT},
) where {FT}
    u_in_tuple = _normalize_velocity(u_in, FT)
    u_sfc_tuple = _normalize_velocity(u_sfc, FT)
    Ts_spec = surface_quantity(Ts, FT)
    qs_spec = surface_quantity(qs, FT)
    z0m_spec = surface_quantity(z0m, FT)
    z0b_spec = surface_quantity(z0b, FT)
    gust_model = _normalize_gustiness(gustiness, FT)
    return SurfaceFluxInputs{
        FT,
        typeof(Ts_spec),
        typeof(qs_spec),
        typeof(z0m_spec),
        typeof(z0b_spec),
        typeof(gust_model),
        typeof(u_in_tuple),
    }(
        Tin,
        qin,
        ρin,
        Ts_spec,
        qs_spec,
        Φs,
        Δz,
        d,
        u_in_tuple,
        u_sfc_tuple,
        gust_model,
        z0m_spec,
        z0b_spec,
        flux_specs.shf,
        flux_specs.lhf,
        flux_specs.ustar,
        flux_specs.Cd,
        flux_specs.Ch,
    )
end

@inline function _normalize_velocity(u::NTuple{2, T}, ::Type{FT}) where {T, FT}
    return (convert(FT, u[1]), convert(FT, u[2]))
end
function _normalize_velocity(u::AbstractVector, ::Type{FT}) where {FT}
    length(u) == 2 ||
        throw(ArgumentError("Velocity vectors must have two horizontal components."))
    return (convert(FT, u[1]), convert(FT, u[2]))
end
_normalize_velocity(u::Nothing, ::Type{FT}) where {FT} = (zero(FT), zero(FT))

@inline function _normalize_gustiness(gustiness, ::Type{FT}) where {FT}
    if gustiness isa GustinessModel
        return gustiness
    elseif gustiness isa Number
        return ConstantGustiness(convert(FT, gustiness))
    else
        return FunctionalGustiness(gustiness)
    end
end

Base.@kwdef mutable struct SurfaceFluxIterationState{FT}
    Ts::FT = FT(0)
    qs::FT = FT(0)
    gustiness::FT = FT(1)
    ustar::FT = FT(0.1)
    L_MO::FT = FT(10)
    shf::FT = FT(0)
    lhf::FT = FT(0)
    Cd::FT = FT(0)
    Ch::FT = FT(0)
    evaporation::FT = FT(0)
    ρ_sfc::FT = FT(1)
    buoyancy_flux::FT = FT(0)
end

struct SimilarityScales{FT}
    u_star::FT
    dsev_star::FT
    q_star::FT
    L_star::FT
    theta_v_star::FT
    ell_u::FT
    ell_theta::FT
    ell_q::FT
end

struct SolverSnapshot{FT, S<:SimilarityScales{FT}}
    scales::S
    ρ_sfc::FT
    gustiness::FT
    Cd::FT
    Ch::FT
    shf::FT
    lhf::FT
    evaporation::FT
    buoyancy_flux::FT
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
    u_in::U
    u_sfc::U
    gustiness::FT
    ustar::FT
    shf::FT
    lhf::FT
    Cd::FT
    Ch::FT
    L_MO::FT
    evaporation::FT
    buoyancy_flux::FT
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
    :u_in,
    :u_sfc,
    :gustiness,
    :ustar,
    :shf,
    :lhf,
    :Cd,
    :Ch,
    :L_MO,
    :evaporation,
    :buoyancy_flux,
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
