"""
    roughness_lengths(momentum; scalar=momentum)

Convenience constructor that pairs the momentum and scalar roughness
specifications.
"""
struct RoughnessLengths{FT}
    momentum::FT
    scalar::FT
end

roughness_lengths(momentum; scalar = momentum) = FixedRoughnessSpec(momentum, scalar)

@inline function default_roughness_lengths(::APS{FT}) where {FT}
    return RoughnessLengths{FT}(FT(1e-3), FT(1e-4))
end

@inline function flux_spec(param_set::APS; kwargs...)
    return FluxSpecs(param_set; kwargs...)
end

gustiness_constant(val::Real) = ConstantGustinessSpec(val)

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
    Ts_guess::FT,
    qs_guess::FT,
    Φs::FT,
    Δz::FT,
    d::FT,
    u_int,
    u_sfc,
    config::SurfaceFluxConfig,
    flux_specs::FluxSpecs{FT},
    update_Ts!,
    update_qs!,
) where {FT}
    roughness_model = instantiate_roughness(param_set, config.roughness)
    gustiness_model = instantiate_gustiness(param_set, config.gustiness)
    return SurfaceFluxInputs(
        Tin,
        qin,
        ρin,
        Ts_guess,
        qs_guess,
        Φs,
        Δz,
        d,
        u_int,
        u_sfc,
        roughness_model,
        gustiness_model,
        update_Ts!,
        update_qs!,
        flux_specs,
    )
end
