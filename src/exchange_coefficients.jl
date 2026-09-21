"""

    drag_coefficient(param_set, ζ, z0m, Δz_eff, scheme, rsl_model = NoRoughnessSubLayer())

Compute the drag coefficient `Cd` for momentum exchange.

# Arguments
- `param_set`: Parameter set
- `ζ`: Stability parameter `ζ = Δz_eff / L_MO`
- `z0m`: Roughness length for momentum [m]
- `Δz_eff`: Effective aerodynamic height `Δz - d` [m]
- `scheme`: Surface flux solver scheme (default: `PointValueScheme()`)
- `rsl_model`: Optional roughness sublayer model (default: [`NoRoughnessSubLayer`](@ref)).

# Formula:

    Cd = (κ / F̂_m)^2

where `F̂_m = F_m + P_m` is the RSL-corrected dimensionless velocity profile
(`P_m ≤ 0` from [`rsl_profile_correction`](@ref)).
"""
function drag_coefficient(
    param_set::APS,
    ζ,
    z0m,
    Δz_eff,
    scheme = UF.PointValueScheme(),
    rsl_model = NoRoughnessSubLayer(),
)
    uf_params = SFP.uf_params(param_set)
    κ = SFP.von_karman_const(param_set)

    F_m =
        UF.dimensionless_profile(uf_params, Δz_eff, ζ, z0m, UF.MomentumTransport(), scheme)
    P_m = rsl_profile_correction(rsl_model, Δz_eff, z0m, UF.MomentumTransport())
    Cd = (κ / (F_m + P_m))^2
    return Cd
end

"""
    drag_coefficient(inputs, speed)

Compute the drag coefficient `Cd` from friction velocity (presumed to be in `inputs.ustar`) 
and effective wind speed (including any gustiness factors).

See the inputs container.

# Arguments
- `inputs`: The inputs container. See [`build_surface_flux_inputs`](@ref SurfaceFluxes.build_surface_flux_inputs).
- `speed`: Effective wind speed [m/s].
"""
function drag_coefficient(inputs, speed)
    ustar = inputs.ustar
    return (ustar / speed)^2
end

"""

    heat_exchange_coefficient(param_set, ζ, z0m, z0h, Δz_eff, scheme, rsl_model = NoRoughnessSubLayer())

Compute the heat exchange coefficient `Ch` for scalar exchange.

# Formula:

    Ch = κ² / (F̂_m · F̂_h),

where `F̂_m = F_m + P_m` and `F̂_h = F_h + P_h` are the RSL-corrected dimensionless
profiles for momentum and scalars respectively (corrections from
[`rsl_profile_correction`](@ref)). For the finite-volume case, this corresponds to
the formulation in Nishizawa & Kitamura (2018), Eqs. 21 & 22 (with Pr_0 absorbed into F_h).

# Arguments
- `param_set`: Parameter set
- `ζ`: Stability parameter `ζ = Δz_eff / L_MO`
- `z0m`: Roughness length for momentum [m]
- `z0h`: Roughness length for scalars (heat/moisture) [m]
- `Δz_eff`: Effective aerodynamic height `Δz - d` [m]
- `scheme`: Surface flux solver scheme (default: `PointValueScheme()`)
- `rsl_model`: Optional roughness sublayer model (default: [`NoRoughnessSubLayer`](@ref)).
"""
function heat_exchange_coefficient(
    param_set::APS,
    ζ,
    z0m,
    z0h,
    Δz_eff,
    scheme = UF.PointValueScheme(),
    rsl_model = NoRoughnessSubLayer(),
)
    uf_params = SFP.uf_params(param_set)
    κ = SFP.von_karman_const(param_set)

    F_m =
        UF.dimensionless_profile(uf_params, Δz_eff, ζ, z0m, UF.MomentumTransport(), scheme)
    F_h = UF.dimensionless_profile(uf_params, Δz_eff, ζ, z0h, UF.HeatTransport(), scheme)
    P_m = rsl_profile_correction(rsl_model, Δz_eff, z0m, UF.MomentumTransport())
    P_h = rsl_profile_correction(rsl_model, Δz_eff, z0h, UF.HeatTransport())

    Ch = κ^2 / ((F_m + P_m) * (F_h + P_h))
    return Ch
end

"""
    heat_conductance(param_set, ζ, ustar, inputs, z0m, z0h, scheme)

Compute the heat conductance `g_h` (speed * Ch), including any gustiness factor in the wind speed.
Calculates windspeed and exchange coefficient internally from Monin-Obukhov variables.


# Arguments
- `param_set`: Parameter set.
- `ζ`: Monin-Obukhov stability parameter.
- `ustar`: Friction velocity [m/s].
- `inputs`: The inputs container. See [`build_surface_flux_inputs`](@ref SurfaceFluxes.build_surface_flux_inputs).
- `z0m`: Momentum roughness length [m].
- `z0h`: Thermal roughness length [m].
- `scheme`: Discretization scheme.
"""
function heat_conductance(
    param_set::APS,
    ζ,
    ustar,
    inputs,
    z0m,
    z0h,
    scheme = UF.PointValueScheme(),
)
    # Compute Ch (pass RSL model from inputs)
    Δz_eff = effective_height(inputs)
    Ch = heat_exchange_coefficient(param_set, ζ, z0m, z0h, Δz_eff, scheme, inputs.rsl_model)

    # Compute windspeed with gustiness (using windspeed helper which handles b_flux)
    current_speed = windspeed(param_set, ζ, ustar, inputs)

    return Ch * current_speed
end
