"""
    build_surface_flux_inputs(args...)

Centralized helper that normalizes user-facing specifications (winds,
roughness, gustiness, flux constraints) into a NamedTuple containing inputs
for surface flux calculations. Input types are preserved (no promotion), and users
should ensure consistent input types for best performance.

# Returns
A `NamedTuple` with the following fields:

## Atmospheric and Surface State
- `T_int`: Interior air temperature [K]
- `q_tot_int`: Interior total specific humidity [kg/kg]
- `q_liq_int`: Interior liquid specific humidity [kg/kg]
- `q_ice_int`: Interior ice specific humidity [kg/kg]
- `ρ_int`: Interior air density [kg/m³]
- `T_sfc_guess`: Initial guess for surface temperature [K]. Can be `nothing` for default fallback.
- `q_vap_sfc_guess`: Initial guess for surface vapor specific humidity [kg/kg]. Can be `nothing` for default fallback.

## Geometry
- `Φ_sfc`: Surface geopotential [m²/s²]
- `Δz`: Height of the reference level above the surface [m], under the convention
  `reference_level`
- `d`: Displacement height [m]
- `reference_level`: Convention for `Δz`, [`ReferenceAboveSurface`](@ref) or
  [`ReferenceAboveApparentSink`](@ref); the solver converts the second to the first
  (see [`reference_above_surface`](@ref))

## Wind
- `u_int`: Horizontal wind components `(u, v)` at the interior level, as a tuple [m/s].
- `u_sfc`: Horizontal wind components `(u, v)` at the surface level, as a tuple [m/s].

## Parameterizations
- `roughness_model`: Roughness parameterization, e.g. [`ConstantRoughnessParams`](@ref).
- `gustiness_model`: Gustiness parameterization, e.g. [`ConstantGustinessSpec`](@ref).
- `moisture_model`: Moisture model, [`MoistModel`](@ref) or [`DryModel`](@ref).
- `rsl_model`: Roughness sublayer model, e.g. [`ExponentialRSL`](@ref) or
  [`NoRoughnessSubLayer`](@ref).
- `stability_cap`: Stability cap specification, e.g. [`MaxHeatFluxStabilityCap`](@ref) or
  [`NoStabilityCap`](@ref).
- `ζ_cap`: Numerical value of the stability cap, or `nothing`. It is `nothing` here and is
  set by the MOST solver from `stability_cap` (see [`with_stability_cap`](@ref)); functions
  that read the cap from the inputs compute it from `stability_cap` when it is `nothing`
  (see [`resolved_stability_cap`](@ref)).
- `roughness_inputs`: Optional inputs for roughness models.

## Callbacks and Prescribed Values
- `update_T_sfc`: Optional callback to update surface temperature during iteration.
- `update_q_vap_sfc`: Optional callback to update surface vapor specific humidity during iteration.
- `shf`: Prescribed sensible heat flux from [`FluxSpecs`](@ref); may be `nothing` [W/m²].
- `lhf`: Prescribed latent heat flux from [`FluxSpecs`](@ref); may be `nothing` [W/m²].
- `ustar`: Prescribed friction velocity from [`FluxSpecs`](@ref); may be `nothing` [m/s].
- `Cd`: Prescribed momentum exchange coefficient from [`FluxSpecs`](@ref); may be `nothing`.
- `Ch`: Prescribed heat exchange coefficient from [`FluxSpecs`](@ref); may be `nothing`.
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
    # Prescribed fluxes and coefficients in the floating-point type of the state, so that
    # Float64 specifications keep a Float32 solve in Float32; `nothing` and dual numbers
    # pass through (see `float_parameter`)
    FT = float(typeof(T_int))

    return (;
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
        u_int = Tuple(u_int),
        u_sfc = Tuple(u_sfc),
        roughness_model = config.roughness,
        gustiness_model = config.gustiness,
        moisture_model = config.moisture_model,
        rsl_model = config.rsl_model,
        stability_cap = config.stability_cap,
        reference_level = config.reference_level,
        ζ_cap = nothing,
        roughness_inputs,
        update_T_sfc,
        update_q_vap_sfc,
        shf = float_parameter(FT, flux_specs.shf),
        lhf = float_parameter(FT, flux_specs.lhf),
        ustar = float_parameter(FT, flux_specs.ustar),
        Cd = float_parameter(FT, flux_specs.Cd),
        Ch = float_parameter(FT, flux_specs.Ch),
    )
end
