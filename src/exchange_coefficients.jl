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
(see [`rsl_corrected_profile`](@ref); `F̂_m = F_m` without RSL).
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

    F̂_m = rsl_corrected_profile(
        uf_params,
        rsl_model,
        Δz_eff,
        ζ,
        z0m,
        UF.MomentumTransport(),
        scheme,
    )
    Cd = (κ / F̂_m)^2
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
profiles for momentum and scalars respectively (see [`rsl_corrected_profile`](@ref)).
For the finite-volume case, this corresponds to the formulation in Nishizawa & Kitamura
(2018), Eqs. 21 & 22 (with Pr_0 absorbed into F_h).

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

    F̂_m = rsl_corrected_profile(
        uf_params,
        rsl_model,
        Δz_eff,
        ζ,
        z0m,
        UF.MomentumTransport(),
        scheme,
    )
    F̂_h =
        rsl_corrected_profile(
            uf_params,
            rsl_model,
            Δz_eff,
            ζ,
            z0h,
            UF.HeatTransport(),
            scheme,
        )

    Ch = κ^2 / (F̂_m * F̂_h)
    return Ch
end

"""
    heat_conductance(param_set, ζ, ustar, inputs, z0m, z0h, scheme)

Compute the heat conductance `g_h = Ch ΔU` [m/s] at stability parameter `ζ`
and friction velocity `ustar`, where `Ch = ϕ_m ϕ_h` is the heat exchange
coefficient at the capped stability parameter ([`MaxHeatFluxStabilityCap`](@ref),
[`capped_stability`](@ref), [`physical_scale_coeff`](@ref)) and `ΔU` is the
effective wind speed including gustiness ([`gustiness_value`](@ref)).

Inside the stability solve, when `inputs` stores the profile coefficients
([`with_profile_coefficients`](@ref)) and the gustiness is independent of
`ustar`, `ustar = ϕ_m ΔU` holds and `heat_conductance` returns `ustar ϕ_h`.

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
    ϕ_h = physical_scale_coeff(param_set, inputs, scheme, ζ, z0h, UF.HeatTransport())
    if ustar_implies_windspeed(inputs)
        return ustar * ϕ_h
    end
    ϕ_m = physical_scale_coeff(param_set, inputs, scheme, ζ, z0m, UF.MomentumTransport())
    current_speed = windspeed(param_set, ζ, ustar, inputs, scheme)
    return ϕ_m * ϕ_h * current_speed
end

# True when `u_star = ϕ_m ΔU` holds with the stored profile coefficient
# `inputs.ϕ_m`: `ϕ_m` is stored, `ustar` is not prescribed, and the gustiness
# model does not depend on `ustar`.
@inline ustar_implies_windspeed(inputs) =
    get(inputs, :ϕ_m, nothing) !== nothing &&
    inputs.ustar === nothing &&
    !depends_on_ustar(inputs.gustiness_model)
