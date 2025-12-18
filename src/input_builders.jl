"""
    gustiness_constant(val)

Create a constant gustiness parameterization with the given value.

Returns a `ConstantGustinessSpec(val)` that provides a fixed gustiness 
contribution to the effective windspeed.
"""
gustiness_constant(val) = ConstantGustinessSpec(val)

# Velocity normalization helper - specialized for 2-tuples to avoid GPU-unsafe throw
@inline _normalize_velocity(u::NTuple{2}, ::Type{FT}) where {FT} = (FT(u[1]), FT(u[2]))
@inline _normalize_velocity(u::AbstractVector, ::Type{FT}) where {FT} = (FT(u[1]), FT(u[2]))


"""
    build_surface_flux_inputs(args...)

Centralized helper that normalizes user-facing specifications (winds,
roughness, gustiness, flux constraints) into a concrete [`SurfaceFluxInputs`](@ref)
instance.
"""
function build_surface_flux_inputs(
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
    # Determine floating point type FT from inputs, handling optional Nothing values
    types = (
        typeof(T_int),
        typeof(q_tot_int),
        typeof(q_liq_int),
        typeof(q_ice_int),
        typeof(ρ_int),
        typeof(Φ_sfc),
        typeof(Δz),
        typeof(d),
    )


    if T_sfc_guess !== nothing
        types = (types..., typeof(T_sfc_guess))
    end
    if q_vap_sfc_guess !== nothing
        types = (types..., typeof(q_vap_sfc_guess))
    end

    FT = promote_type(types...)

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
        typeof(flux_specs.shf),
        typeof(flux_specs.lhf),
        typeof(flux_specs.ustar),
        typeof(flux_specs.Cd),
        typeof(flux_specs.Ch),
    }(
        FT(T_int),
        FT(q_tot_int),
        FT(q_liq_int),
        FT(q_ice_int),
        FT(ρ_int),
        T_sfc_guess === nothing ? nothing : FT(T_sfc_guess),
        q_vap_sfc_guess === nothing ? nothing : FT(q_vap_sfc_guess),
        FT(Φ_sfc),
        FT(Δz),
        FT(d),
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
