"""
    gustiness_value(spec, param_set, buoyancy_flux)

Returns the gustiness velocity scale [m/s] based on the specification.

# Arguments
- `spec`: The gustiness specification (e.g., [`ConstantGustinessSpec`](@ref) or [`DeardorffGustinessSpec`](@ref)).
- `param_set`: Parameter set containing constants and coefficients.
- `buoyancy_flux`: Surface buoyancy flux [m^2/s^3], required for Deardorff gustiness.

"""
@inline gustiness_value(spec::ConstantGustinessSpec, param_set, buoyancy_flux) =
    float_parameter(eltype(param_set), spec.value)

"""
    gustiness_value(::DeardorffGustinessSpec, param_set, buoyancy_flux)

Calculates the gustiness based on the Deardorff convective velocity scale.

# Formulation
The gustiness ``U_{gust}`` is parameterized as proportional to the Deardorff velocity ``w_*``:
```math
U_{gust} = C_{gust} \\cdot w_*
```
where ``w_* = (B \\cdot z_i)^{1/3}``.

- ``B`` is the surface buoyancy flux (`buoyancy_flux`).
- ``z_i`` is the boundary layer height (assumed fixed to`gustiness_zi` from parameters).
- ``C_{gust}`` is a scaling coefficient (`gustiness_coeff` from parameters).

This formulation parametrizes the enhancement of surface fluxes due to boundary layer scale
eddies in unstable conditions, particularly important in low-wind regimes
(free convection limit).

# References
- Deardorff, J. W. (1970). Convective velocity and temperature scales for the unstable planetary
  boundary layer and for Rayleigh convection. Journal of the Atmospheric Sciences, 27, 1211-1213.
  [DOI: 10.1175/1520-0469(1970)027<1211:CVATSF>2.0.CO;2](https://doi.org/10.1175/1520-0469(1970)027%3C1211:CVATSF%3E2.0.CO;2)
- Beljaars, A. C. M. (1995). The parametrization of surface fluxes in large-scale models under free convection 
  Quarterly Journal of the Royal Meteorological Society, 121, 255-270.
  [DOI:  10.1002/qj.49712152203](https://doi.org/10.1002/qj.49712152203)
"""
@inline function gustiness_value(::DeardorffGustinessSpec, param_set, buoyancy_flux)
    FT = eltype(param_set)
    # Extract parameters
    β = SFP.gustiness_coeff(param_set)
    zi = SFP.gustiness_zi(param_set)

    w_star = cbrt(max(buoyancy_flux * zi, FT(0)))
    return β * w_star
end

"""
    Δu_components(inputs)

Computes the vector difference between the interior and surface wind components.

Returns a tuple `(Δu_x, Δu_y)`.
"""
@inline function Δu_components(inputs)
    return (
        inputs.u_int[1] - inputs.u_sfc[1],
        inputs.u_int[2] - inputs.u_sfc[2],
    )
end

"""
    windspeed(Δu, gustiness)
    windspeed(inputs, gustiness)

Computes the effective wind speed magnitude [m/s], accounting for gustiness.

# Formulation
The effective wind speed is calculated as the maximum of the mean wind speed difference
and the gustiness scale:
```math
U_{\\text{eff}} = \\max(\\sqrt{\\Delta u_x^2 + \\Delta u_y^2}, U_{gust})
```
This formulation ensures that surface fluxes remain non-zero even in the absence of mean wind,
driven by convective eddies or other sub-grid variability represented by ``U_{gust}``. This is 
important in low-wind regimes in the free convection limit.

# Arguments
- `Δu`: Tuple of wind component differences `(Δu_x, Δu_y)` [m/s]. The components may
  differ in type, e.g., a dual number and a float when differentiating with respect to
  one component.
- `inputs`: The inputs container. See [`build_surface_flux_inputs`](@ref SurfaceFluxes.build_surface_flux_inputs).
- `gustiness`: Gustiness velocity scale [m/s].
"""
@inline function windspeed(Δu::Tuple{Any, Any}, gustiness)
    return max(hypot(Δu[1], Δu[2]), gustiness)
end

@inline function windspeed(inputs, gustiness)
    return windspeed(Δu_components(inputs), gustiness)
end

"""
    windspeed(inputs, param_set, buoyancy_flux)

Computes the effective wind speed magnitude [m/s], including any gustiness factor.
"""
@inline function windspeed(inputs, param_set, buoyancy_flux)
    gustiness = gustiness_value(inputs.gustiness_model, param_set, buoyancy_flux)
    return windspeed(inputs, gustiness)
end

"""
    gustiness_value(spec, param_set, ζ, ustar, inputs, scheme = PointValueScheme())

Return the gustiness velocity scale [m/s] from the solver variables `ζ` and `ustar`.
[`ConstantGustinessSpec`](@ref) returns its value; [`FlooredDeardorffGustinessSpec`](@ref)
evaluates its convective part in closed form with the profile integrals of `scheme`; the
other models evaluate the buoyancy flux first (see [`depends_on_ustar`](@ref)).
"""
@inline gustiness_value(
    spec::ConstantGustinessSpec,
    param_set,
    ζ,
    ustar,
    inputs,
    scheme = PointValueScheme(),
) = float_parameter(eltype(param_set), spec.value)
@inline function gustiness_value(
    spec::AbstractGustinessSpec,
    param_set,
    ζ,
    ustar,
    inputs,
    scheme = PointValueScheme(),
)
    b_flux = buoyancy_flux(param_set, ζ, ustar, inputs)
    return gustiness_value(spec, param_set, b_flux)
end

"""
    gustiness_value(spec::FlooredDeardorffGustinessSpec, param_set, buoyancy_flux)

Return the larger of the floor `spec.u_min` and the convective gustiness ``β w_*`` of
[`DeardorffGustinessSpec`](@ref) for the surface buoyancy flux `buoyancy_flux` [m²/s³],
zero for a non-positive buoyancy flux. The post-solve fluxes use this form, with the
buoyancy flux implied by the converged `ζ` and friction velocity.
"""
@inline gustiness_value(spec::FlooredDeardorffGustinessSpec, param_set, buoyancy_flux) =
    max(
        float_parameter(eltype(param_set), spec.u_min),
        gustiness_value(DeardorffGustinessSpec(), param_set, buoyancy_flux),
    )

"""
    gustiness_value(
        spec::FlooredDeardorffGustinessSpec, param_set, ζ, ustar, inputs,
        scheme = PointValueScheme(),
    )
    gustiness_value(spec, param_set, ζ, ustar, inputs, scheme, Δθ_v, θ_v_int)

Return the larger of the floor `spec.u_min` and the free-convection wind speed
[`free_convection_wind_speed`](@ref) at the stability parameter `ζ`, with the profile
integrals of `scheme`. The friction velocity `ustar` enters only through the roughness
lengths. The second form takes the virtual potential temperature excess `Δθ_v` of the
surface over the interior air and the interior value `θ_v_int` [K] from the caller; for
the other gustiness models, it ignores them.
"""
@inline function gustiness_value(
    spec::FlooredDeardorffGustinessSpec,
    param_set,
    ζ,
    ustar,
    inputs,
    scheme = PointValueScheme(),
)
    return max(
        float_parameter(eltype(param_set), spec.u_min),
        free_convection_wind_speed(param_set, ζ, ustar, inputs, scheme),
    )
end
@inline function gustiness_value(
    spec::FlooredDeardorffGustinessSpec,
    param_set,
    ζ,
    ustar,
    inputs,
    scheme,
    Δθ_v,
    θ_v_int,
)
    return max(
        float_parameter(eltype(param_set), spec.u_min),
        free_convection_wind_speed(param_set, ζ, ustar, inputs, scheme, Δθ_v, θ_v_int),
    )
end
@inline gustiness_value(
    spec::AbstractGustinessSpec,
    param_set,
    ζ,
    ustar,
    inputs,
    scheme,
    Δθ_v,
    θ_v_int,
) = gustiness_value(spec, param_set, ζ, ustar, inputs, scheme)

"""
    free_convection_wind_speed(param_set, ζ, ustar, inputs, scheme = PointValueScheme())
    free_convection_wind_speed(param_set, ζ, ustar, inputs, scheme, Δθ_v, θ_v_int)

Return the effective wind speed ``U`` [m/s] at which the convective gustiness ``β w_*``
equals ``U`` itself, for the Monin-Obukhov stability parameter `ζ` and the surface and
atmospheric state in `inputs`:

```math
U^2 = β^3 κ^2 \\frac{g}{θ_v} z_i \\frac{Δθ_v}{F_m(ζ) F_h(ζ)},
```

where ``Δθ_v`` is the virtual potential temperature excess of the surface over the air,
``θ_v`` the air value, ``F_m`` and ``F_h`` the dimensionless profile integrals of momentum
and heat (so that ``u_* = κ U / F_m`` and ``θ_{v*} = κ Δθ_v / F_h``), and ``β``, ``z_i``,
``κ``, ``g`` the parameters `gustiness_coeff`, `gustiness_zi`, `von_karman_const`, and
`grav`. The result follows from ``w_*^3 = B z_i`` with ``B = (g / θ_v) u_* θ_{v*}`` and
``U = β w_*``. It is zero in stable conditions (``ζ > 0``) or when ``Δθ_v ≤ 0``,
matching the post-solve gustiness from the buoyancy flux ``-u_*^3 ζ / (κ Δz_{eff})``.
The profile integrals are evaluated with the discretization `scheme` of the solve
([`PointValueScheme`](@ref) or [`LayerAverageScheme`](@ref)), the roughness sublayer
model, and the stability cap of `inputs`, with the roughness lengths of the roughness
model at `ustar`. The surface temperature and humidity are `T_sfc_guess` and
`q_vap_sfc_guess` of `inputs` (updated by the surface callbacks when present), so the
gustiness and the state bulk Richardson number see the same surface state; the second
form takes `Δθ_v` and `θ_v_int` [K] ([`virtual_pottemps`](@ref)) from the caller.
"""
@inline function free_convection_wind_speed(
    param_set::APS,
    ζ,
    ustar,
    inputs,
    scheme = PointValueScheme(),
)
    T_sfc = safe_T_sfc_guess(inputs)
    q_vap_sfc = safe_q_vap_sfc_guess(inputs)
    ρ_sfc = surface_density(param_set, inputs, T_sfc, q_vap_sfc)
    θ_v_sfc, θ_v_int = virtual_pottemps(param_set, inputs, T_sfc, ρ_sfc, q_vap_sfc)
    Δθ_v = θ_v_sfc - θ_v_int
    return free_convection_wind_speed(param_set, ζ, ustar, inputs, scheme, Δθ_v, θ_v_int)
end
@inline function free_convection_wind_speed(
    param_set::APS,
    ζ,
    ustar,
    inputs,
    scheme,
    Δθ_v,
    θ_v_int,
)
    β = SFP.gustiness_coeff(param_set)
    z_i = SFP.gustiness_zi(param_set)
    g = SFP.grav(param_set)
    z0m, z0h = roughness_lengths(param_set, inputs, ustar)
    # κ / F_m and κ / F_h
    ϕ_m = physical_scale_coeff(param_set, inputs, scheme, ζ, z0m, UF.MomentumTransport())
    ϕ_h = physical_scale_coeff(param_set, inputs, scheme, ζ, z0h, UF.HeatTransport())
    U²_raw = β^3 * (g / θ_v_int) * z_i * Δθ_v * ϕ_m * ϕ_h
    U² = ifelse(ζ > zero(ζ), zero(U²_raw), U²_raw)
    return sqrt(max(U², zero(U²)))
end

"""
    windspeed(param_set, ζ, ustar, inputs, scheme = PointValueScheme())

Compute the effective wind speed magnitude [m/s] from the solver variables `ζ` and
`ustar`, with the gustiness from [`gustiness_value`](@ref) at the profile integrals of
`scheme`.
"""
@inline function windspeed(param_set::APS, ζ, ustar, inputs, scheme = PointValueScheme())
    gustiness =
        gustiness_value(inputs.gustiness_model, param_set, ζ, ustar, inputs, scheme)
    return windspeed(inputs, gustiness)
end

"""
    minimum_wind_speed(spec::AbstractGustinessSpec, param_set)

Return the minimum effective wind speed [m/s] that a gustiness model imposes in all
conditions, in the floating-point type of `param_set`: the value of a
[`ConstantGustinessSpec`](@ref), the floor `u_min` of a
[`FlooredDeardorffGustinessSpec`](@ref), and zero for models whose gustiness vanishes in
stable conditions, such as [`DeardorffGustinessSpec`](@ref). A model that folds this
floor into the wind speed it passes to the solve (for example, a canopy model that
attenuates the wind above the canopy to the ground below it) pairs the attenuated wind
with [`without_floor`](@ref).
"""
@inline minimum_wind_speed(spec::ConstantGustinessSpec, param_set::APS) =
    float_parameter(eltype(param_set), spec.value)
@inline minimum_wind_speed(spec::FlooredDeardorffGustinessSpec, param_set::APS) =
    float_parameter(eltype(param_set), spec.u_min)
@inline minimum_wind_speed(::AbstractGustinessSpec, param_set::APS) =
    zero(eltype(param_set))

"""
    without_floor(spec::AbstractGustinessSpec)

Return the gustiness model `spec` with its minimum wind speed set to zero: a
[`ConstantGustinessSpec`](@ref) becomes a zero gustiness, a
[`FlooredDeardorffGustinessSpec`](@ref) keeps its convective part only, and models with
no floor are returned as they are. See [`minimum_wind_speed`](@ref).
"""
@inline without_floor(spec::ConstantGustinessSpec) = ConstantGustinessSpec(zero(spec.value))
@inline without_floor(spec::FlooredDeardorffGustinessSpec) =
    FlooredDeardorffGustinessSpec(zero(spec.u_min))
@inline without_floor(spec::AbstractGustinessSpec) = spec

"""
    wind_depends_on_surface_state(spec::AbstractGustinessSpec)

Whether the effective wind speed of a gustiness model at fixed stability parameter `ζ`
depends on the surface temperature and humidity. The gustiness of
[`FlooredDeardorffGustinessSpec`](@ref) depends on the virtual potential temperature
excess of the surface; that of [`ConstantGustinessSpec`](@ref) is fixed, and that of
[`DeardorffGustinessSpec`](@ref) follows from the buoyancy flux implied by `ζ` and the
friction velocity. For models independent of the surface state, the surface callbacks
take a single step at each `ζ`, and the wind speed, friction velocity, and roughness
lengths are not re-evaluated after the step.
"""
@inline wind_depends_on_surface_state(::AbstractGustinessSpec) = true
@inline wind_depends_on_surface_state(::ConstantGustinessSpec) = false
@inline wind_depends_on_surface_state(::DeardorffGustinessSpec) = false
