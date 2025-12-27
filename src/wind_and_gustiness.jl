"""
    gustiness_value(spec, param_set, buoyancy_flux)

Returns the gustiness velocity scale [m/s] based on the specification.

# Arguments
- `spec`: The gustiness specification (e.g., [`ConstantGustinessSpec`](@ref) or [`DeardorffGustinessSpec`](@ref)).
- `param_set`: Parameter set containing constants and coefficients.
- `buoyancy_flux`: Surface buoyancy flux [m^2/s^3], required for Deardorff gustiness.

"""
@inline gustiness_value(spec::ConstantGustinessSpec, param_set, buoyancy_flux) = spec.value

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
    # Extract parameters
    β = SFP.gustiness_coeff(param_set)
    zi = SFP.gustiness_zi(param_set)

    w_star = cbrt(max(buoyancy_flux * zi, 0))
    return β * w_star
end

"""
    Δu_components(inputs::SurfaceFluxInputs)

Computes the vector difference between the interior and surface wind components.

Returns a tuple `(Δu_x, Δu_y)`.
"""
@inline function Δu_components(inputs::SurfaceFluxInputs)
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
U_{\text{eff}} = \\max(\\sqrt{\\Delta u_x^2 + \\Delta u_y^2}, U_{gust})
```
This formulation ensures that surface fluxes remain non-zero even in the absence of mean wind,
driven by convective eddies or other sub-grid variability represented by ``U_{gust}``. This is 
important in low-wind regimes in the free convection limit.

# Arguments
- `Δu`: Tuple of wind component differences `(Δu_x, Δu_y)`.
- `inputs`: [`SurfaceFluxInputs`](@ref) struct (convenience wrapper).
- `gustiness`: Gustiness velocity scale [m/s].
"""
@inline function windspeed(Δu::NTuple{2}, gustiness)
    return max(hypot(Δu[1], Δu[2]), gustiness)
end

@inline function windspeed(inputs::SurfaceFluxInputs, gustiness)
    return windspeed(Δu_components(inputs), gustiness)
end

"""
    windspeed(inputs, param_set, buoyancy_flux)

Computes the effective wind speed magnitude [m/s], including any gustiness factor.
"""
@inline function windspeed(inputs::SurfaceFluxInputs, param_set, buoyancy_flux)
    gustiness = gustiness_value(inputs.gustiness_model, param_set, buoyancy_flux)
    return windspeed(inputs, gustiness)
end

"""
    windspeed(param_set, ζ, ustar, inputs)

Computes the effective wind speed magnitude [m/s] from solver variables.
Calculates buoyancy flux and gustiness internally from Monin-Obukhov variables.
"""
@inline function windspeed(param_set::APS, ζ, ustar, inputs::SurfaceFluxInputs)
    b_flux = buoyancy_flux(param_set, ζ, ustar, inputs)
    gustiness = gustiness_value(inputs.gustiness_model, param_set, b_flux)
    return windspeed(inputs, gustiness)
end
