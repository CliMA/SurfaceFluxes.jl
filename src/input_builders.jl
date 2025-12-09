

@inline function flux_spec(param_set::APS; kwargs...)
    return FluxSpecs(param_set; kwargs...)
end

gustiness_constant(val) = ConstantGustinessSpec(val)

# Velocity normalization helper (moved from types.jl)
@inline function _normalize_velocity(u::Tuple{Any, Any}, ::Type{FT}) where {FT}
    return (convert(FT, u[1]), convert(FT, u[2]))
end
function _normalize_velocity(u::AbstractVector, ::Type{FT}) where {FT}
    length(u) == 2 ||
        throw(ArgumentError("Velocity vectors must have two horizontal components."))
    return (convert(FT, u[1]), convert(FT, u[2]))
end
_normalize_velocity(u::Nothing, ::Type{FT}) where {FT} = (zero(FT), zero(FT))
_normalize_velocity(u::Tuple{}, ::Type{FT}) where {FT} = (zero(FT), zero(FT))


"""
    build_surface_flux_inputs(param_set, args...)

Centralized helper that normalizes user-facing specifications (winds,
roughness, gustiness, flux constraints) into a concrete `SurfaceFluxInputs`
instance.
"""
function build_surface_flux_inputs(
    param_set::APS{FT},
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
    config::SurfaceFluxConfig,
    roughness_inputs,
    flux_specs::FluxSpecs{FT},
    update_Ts!,
    update_qs!,
) where {FT}
    u_int_tuple = _normalize_velocity(u_int, FT)
    u_sfc_tuple = _normalize_velocity(u_sfc, FT)
    
    return SurfaceFluxInputs{
        FT,
        typeof(config.roughness),
        typeof(config.gustiness),
        typeof(roughness_inputs),
        typeof(update_Ts!),
        typeof(update_qs!),
        typeof(u_int_tuple)
    }(
        convert(FT, Tin),
        convert(FT, qin),
        convert(FT, ρin),
        convert(FT, Ts_guess),
        convert(FT, qs_guess),
        convert(FT, Φs),
        convert(FT, Δz),
        convert(FT, d),
        u_int_tuple,
        u_sfc_tuple,
        config.roughness,
        config.gustiness,
        roughness_inputs,
        update_Ts!,
        update_qs!,
        flux_specs.shf,
        flux_specs.lhf,
        flux_specs.ustar,
        flux_specs.Cd,
        flux_specs.Ch,
    )
end
