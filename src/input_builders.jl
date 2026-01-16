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
- `Δz`: Height difference between interior and surface reference levels [m]
- `d`: Displacement height [m]

## Wind
- `u_int`: Horizontal wind components (u, v) at interior level [m/s] (tuple)
- `u_sfc`: Horizontal wind components (u, v) at surface level [m/s] (tuple)

## Parameterizations
- `roughness_model`: Roughness parameterization (e.g., [`ConstantRoughnessParams`](@ref))
- `gustiness_model`: Gustiness parameterization (e.g., [`ConstantGustinessSpec`](@ref))
- `moisture_model`: Moisture model ([`MoistModel`](@ref) or [`DryModel`](@ref))
- `roughness_inputs`: Optional inputs for roughness models

## Callbacks and Prescribed Values
- `update_T_sfc`: Optional callback to update surface temperature during iteration
- `update_q_vap_sfc`: Optional callback to update surface vapor specific humidity during iteration
- `shf`: Prescribed sensible heat flux [W/m²] (from [`FluxSpecs`](@ref), can be `nothing`)
- `lhf`: Prescribed latent heat flux [W/m²] (from [`FluxSpecs`](@ref), can be `nothing`)
- `ustar`: Prescribed friction velocity [m/s] (from [`FluxSpecs`](@ref), can be `nothing`)
- `Cd`: Prescribed momentum exchange coefficient (from [`FluxSpecs`](@ref), can be `nothing`)
- `Ch`: Prescribed heat exchange coefficient (from [`FluxSpecs`](@ref), can be `nothing`)
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
        roughness_inputs,
        update_T_sfc,
        update_q_vap_sfc,
        shf = flux_specs.shf,
        lhf = flux_specs.lhf,
        ustar = flux_specs.ustar,
        Cd = flux_specs.Cd,
        Ch = flux_specs.Ch,
    )
end
