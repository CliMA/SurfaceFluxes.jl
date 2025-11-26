"""
    roughness_lengths(momentum; scalar=momentum)

Convenience constructor that pairs the momentum and scalar roughness
specifications. `momentum` and `scalar` may be scalars, callables, or any of
the `SurfaceQuantity` helper specs.
"""
struct RoughnessLengths{ZM, ZB}
    momentum::ZM
    scalar::ZB
end

roughness_lengths(momentum; scalar = momentum) = RoughnessLengths(momentum, scalar)

@inline function default_roughness_lengths(::APS{FT}) where {FT}
    return RoughnessLengths(FT(1e-3), FT(1e-4))
end

@inline function normalize_roughness(roughness::RoughnessLengths, ::Type{FT}) where {FT}
    return roughness
end

@inline function normalize_roughness(
    roughness::NTuple{2, Any},
    ::Type{FT},
) where {FT}
    return RoughnessLengths(roughness[1], roughness[2])
end

@inline function normalize_roughness(roughness::NamedTuple, ::Type{FT}) where {FT}
    first_value = first(values(roughness))
    momentum = get(roughness, :momentum, get(roughness, :scalar, first_value))
    scalar = get(roughness, :scalar, momentum)
    return RoughnessLengths(momentum, scalar)
end

@inline function normalize_roughness(roughness, ::Type{FT}) where {FT}
    return RoughnessLengths(roughness, roughness)
end

@inline function flux_spec(param_set::APS; kwargs...)
    return FluxSpecs(param_set; kwargs...)
end

gustiness_constant(val::Real) = ConstantGustiness(val)

@inline gustiness_callable(fn) = FunctionalGustiness(fn)

"""
    build_surface_flux_inputs(param_set, args...)

Centralized helper that normalizes user-facing specifications (winds,
roughness, gustiness, flux constraints) into a concrete `SurfaceFluxInputs`
instance. This keeps the GPU-facing solver free of keyword arguments while
presents a clean API to callers.
"""
function build_surface_flux_inputs(
    param_set::APS{FT},
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
    roughness,
    gustiness,
    flux_specs::FluxSpecs{FT},
) where {FT}
    normalized_roughness = normalize_roughness(roughness, FT)
    return SurfaceFluxInputs(
        Tin,
        qin,
        ρin,
        Ts,
        qs,
        Φs,
        Δz,
        d,
        u_in,
        u_sfc,
        gustiness,
        normalized_roughness.momentum,
        normalized_roughness.scalar,
        flux_specs,
    )
end

