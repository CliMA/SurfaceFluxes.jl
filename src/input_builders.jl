

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
    param_set::APS,
    T_int,
    q_tot_int,
    q_liq_int,
    q_ice_int,
    ρ_int,
    T_sfc_guess,
    q_vap_sfc_guess,
    Φ_sfc,
    Δz,
    d,
    u_int,
    u_sfc,
    config::SurfaceFluxConfig,
    roughness_inputs,
    flux_specs,
    update_T_sfc,
    update_q_vap_sfc,
)
    FT = promote_type(
        typeof(T_int),
        typeof(q_tot_int),
        typeof(q_liq_int),
        typeof(q_ice_int),
        typeof(ρ_int),
        typeof(T_sfc_guess),
        typeof(q_vap_sfc_guess),
        typeof(Φ_sfc),
        typeof(Δz),
        typeof(d)
    )
    u_int_tuple = _normalize_velocity(u_int, FT)
    u_sfc_tuple = _normalize_velocity(u_sfc, FT)
    
    return SurfaceFluxInputs{
        FT,
        typeof(config.roughness),
        typeof(config.gustiness),
        typeof(config.moisture_model),
        typeof(roughness_inputs),
        typeof(update_T_sfc),
        typeof(update_q_vap_sfc),
    }(
        convert(FT, T_int),
        convert(FT, q_tot_int),
        convert(FT, q_liq_int),
        convert(FT, q_ice_int),
        convert(FT, ρ_int),
        convert(FT, T_sfc_guess),
        convert(FT, q_vap_sfc_guess),
        convert(FT, Φ_sfc),
        convert(FT, Δz),
        convert(FT, d),
        u_int_tuple,
        u_sfc_tuple,
        config.roughness,
        config.gustiness,
        config.moisture_model,
        roughness_inputs,
        update_T_sfc,
        update_q_vap_sfc,
        flux_specs.shf,
        flux_specs.lhf,
        flux_specs.ustar,
        flux_specs.Cd,
        flux_specs.Ch,
    )
end
