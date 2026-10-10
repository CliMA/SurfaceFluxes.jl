"""
    compute_profile_value(param_set, L_MO, z0, Δz_eff, scale, val_sfc, transport, scheme, rsl_model)

Compute the (nondimensional) value of a variable (momentum or scalar) 
at effective aerodynamic height `Δz_eff` (height above surface minus displacement height).

# Arguments
- `param_set`: Parameter set.
- `L_MO`: Monin-Obukhov length [m].
- `z0`: Roughness length [m].
- `Δz_eff`: Effective aerodynamic height `z - d` [m].
- `scale`: Similarity scale (u_star, theta_star, etc.).
- `val_sfc`: Surface value of the variable.
- `transport`: Transport type (`MomentumTransport` or `HeatTransport`).
- `scheme`: Discretization scheme (default: `PointValueScheme()`).
- `rsl_model`: Roughness sublayer model (default: [`NoRoughnessSubLayer`](@ref)). Use the
  same model as in the flux calculation for consistent profiles.

# Formula:

    X(Δz_eff) = (scale / κ) * F̂_z + val_sfc

where `F̂_z = F_z + P` is the dimensionless profile at height `Δz_eff`, including the
roughness sublayer correction `P` (see [`rsl_corrected_profile`](@ref)).

!!! warning "Stability caps"
    With a stability cap (e.g., [`MaxHeatFluxStabilityCap`](@ref)), the returned `L_MO`
    is the Obukhov length implied by the fluxes, but the exchange coefficients and
    similarity scales were evaluated at the capped stability parameter
    `min(ζ, ζ_cap)` at the forcing height `Δz_eff_ref`. Profiles consistent with the
    fluxes (which reproduce the forcing values at `Δz_eff_ref`) are obtained by passing
    the effective length `L_eff = Δz_eff_ref / min(ζ, ζ_cap)`, returned as the field
    `L_eff` of [`SurfaceFluxConditions`](@ref), instead of `L_MO`. Passing `L_MO` beyond
    the cap overestimates the recovered differences.
"""
function compute_profile_value(
    param_set::APS,
    L_MO,
    z0,
    Δz_eff,
    scale,
    val_sfc,
    transport,
    scheme = UF.PointValueScheme(),
    rsl_model = NoRoughnessSubLayer(),
)
    uf_params = SFP.uf_params(param_set)
    κ = SFP.von_karman_const(param_set)
    ζ = Δz_eff / L_MO

    F̂ = rsl_corrected_profile(uf_params, rsl_model, Δz_eff, ζ, z0, transport, scheme)

    return F̂ * scale / κ + val_sfc
end

"""
    dimensionless_profile_value(param_set, L_eff, z0, z, Δz_eff, transport, scheme, rsl_model)

Return the dimensionless Monin-Obukhov profile ``\\widehat{F}(z)`` of
[`compute_profile_value`](@ref) at the height `z` above the displacement height, for the
Obukhov length `L_eff` of a solve at the reference height `Δz_eff`. The height is
clamped to the range from `z0`, where the profile is zero, to `Δz_eff`, so that the
profile stays within the levels the solve connected.

# Arguments
- `param_set`: Parameter set.
- `L_eff`: Effective Obukhov length of the solve, the field `L_eff` of
  [`SurfaceFluxConditions`](@ref) [m].
- `z0`: Roughness length of the transported quantity [m].
- `z`: Height above the displacement height [m].
- `Δz_eff`: Height of the reference level above the displacement height [m].
- `transport`: `UF.MomentumTransport()` or `UF.HeatTransport()`.
- `scheme`: Discretization scheme ([`PointValueScheme`](@ref) or
  [`LayerAverageScheme`](@ref)).
- `rsl_model`: Roughness sublayer model of the solve (e.g., [`NoRoughnessSubLayer`](@ref)).
"""
@inline function dimensionless_profile_value(
    param_set::APS,
    L_eff,
    z0,
    z,
    Δz_eff,
    transport,
    scheme,
    rsl_model,
)
    FT = eltype(param_set)
    κ = SFP.von_karman_const(param_set)
    z_clamped = max(min(float_parameter(FT, z), Δz_eff), z0)
    # With scale κ and zero surface value, the profile value is F̂ itself
    return compute_profile_value(
        param_set,
        L_eff,
        z0,
        z_clamped,
        κ,
        zero(FT),
        transport,
        scheme,
        rsl_model,
    )
end

"""
    screen_level_values(param_set, sc, inputs, z_screen, z_anemometer, scheme = PointValueScheme())

Return the NamedTuple `(; T, q, u)` of the air temperature [K] and vapor specific
humidity [kg/kg] at the height `z_screen` [m] above the apparent sink for heat
`d + z0h`, and the wind speed [m/s] at the height `z_anemometer` [m] above the apparent
sink for momentum `d + z0m`, reconstructed from the Monin-Obukhov profiles of the solve
that returned `sc` for `inputs`. The WMO screen and anemometer heights are 2 m and 10 m.

Between the surface and the reference level, a scalar that follows the heat profile takes
the value ``X_{sfc} + (X_{int} - X_{sfc}) r(z)`` with
``r(z) = \\widehat{F}_h(z) / \\widehat{F}_h(Δz_{eff})`` (see
[`dimensionless_profile_value`](@ref)), evaluated at the effective Obukhov length
`sc.L_eff`, so that the profile reproduces the fluxes also under a stability cap. The
temperature follows this relation in terms of the dry static energy, the variable the
sensible heat flux is computed from, with the surface state at the displacement height
(see [`surface_geopotential`](@ref)), and so includes the adiabatic change
``g / c_{p,d}`` per meter between the screen and reference levels. The wind speed is
``u_* \\widehat{F}_m(z) / κ``, gustiness included, relative to the surface velocity
`u_sfc` of the inputs.

The screen and anemometer values are point values of the profiles, while the reference
profile ``\\widehat{F}_h(Δz_{eff})`` follows `scheme`. Under [`PointValueScheme`](@ref),
the profiles reach the interior state at the reference level, and the wind reaches the
effective wind speed of the solve. Under [`LayerAverageScheme`](@ref), the interior state is
a layer average, which the point profile attains low in the layer, so point values in
the upper part of the layer lie farther from the surface values than the interior state.
Levels above the reference level take the values at the reference level. At or below
the roughness length, the apparent sink at `d + z0` above the surface, the dry static
energy and humidity take the surface values and the wind vanishes, unless a
roughness-sublayer model is configured: its correction keeps the profiles away from the
surface values there.

# Arguments
- `param_set`: Parameter set.
- `sc`: The [`SurfaceFluxConditions`](@ref) returned by [`surface_fluxes`](@ref).
- `inputs`: The inputs container of that solve. See [`build_surface_flux_inputs`](@ref SurfaceFluxes.build_surface_flux_inputs).
- `z_screen`: Height of the screen level above the apparent sink for heat [m].
- `z_anemometer`: Height of the anemometer above the apparent sink for momentum [m].
- `scheme`: Discretization scheme of the solve (default: [`PointValueScheme`](@ref)).

# Returns
A NamedTuple `(; T, q, u)`:
- `T`: Air temperature at the screen height [K].
- `q`: Vapor specific humidity at the screen height [kg/kg].
- `u`: Wind speed at the anemometer height, relative to `u_sfc` [m/s].
"""
@inline function screen_level_values(
    param_set::APS,
    sc::SurfaceFluxConditions,
    inputs,
    z_screen,
    z_anemometer,
    scheme = UF.PointValueScheme(),
)
    FT = eltype(param_set)
    thermo_params = SFP.thermodynamics_params(param_set)
    κ = SFP.von_karman_const(param_set)
    g = SFP.grav(param_set)
    cp_d = TD.Parameters.cp_d(thermo_params)
    inputs = reference_above_surface(param_set, inputs)
    Δz_eff = effective_height(param_set, inputs)
    z0m, z0h = momentum_and_scalar_roughness(
        inputs.roughness_model,
        sc.ustar,
        param_set,
        inputs.roughness_inputs,
    )
    heat = UF.HeatTransport()
    momentum = UF.MomentumTransport()
    rsl = inputs.rsl_model

    # Heights above the displacement height of the screen level and of the anemometer
    z_T = z0h + float_parameter(FT, z_screen)
    z_u = z0m + float_parameter(FT, z_anemometer)

    # The reference profile follows the scheme of the solve, which connects the interior
    # values (layer averages under LayerAverageScheme) to the surface; the screen and
    # anemometer values are point values
    point = UF.PointValueScheme()
    F̂_h_ref = dimensionless_profile_value(
        param_set,
        sc.L_eff,
        z0h,
        Δz_eff,
        Δz_eff,
        heat,
        scheme,
        rsl,
    )
    F̂_h = dimensionless_profile_value(
        param_set,
        sc.L_eff,
        z0h,
        z_T,
        Δz_eff,
        heat,
        point,
        rsl,
    )
    r = ifelse(F̂_h_ref > 0, F̂_h / max(F̂_h_ref, eps(FT)), FT(1))

    # The dry static energy varies linearly with r between the surface state, at the
    # displacement height, and the reference level; the heights above the displacement
    # height of the reference and screen levels convert it back to temperature
    T_sfc = sc.T_sfc
    T_int = inputs.T_int
    z_T_clamped = max(min(z_T, Δz_eff), z0h)
    T = T_sfc + (T_int - T_sfc) * r + g / cp_d * (r * Δz_eff - z_T_clamped)
    q_sfc = sc.q_vap_sfc
    q = q_sfc + (interior_vapor_specific_humidity(inputs) - q_sfc) * r

    F̂_m = dimensionless_profile_value(
        param_set,
        sc.L_eff,
        z0m,
        z_u,
        Δz_eff,
        momentum,
        point,
        rsl,
    )
    u = sc.ustar * max(F̂_m, FT(0)) / κ
    return (; T, q, u)
end
