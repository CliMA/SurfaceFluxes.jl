"""
    drag_coefficient(param_set, ζ, z0m, Δz, scheme)

Compute the drag coefficient `Cd` for momentum exchange.

# Arguments
- `param_set`: Parameter set
- `ζ`: Stability parameter `ζ = Δz / L_MO`
- `z0m`: Roughness length for momentum [m]
- `Δz`: Height difference between the surface and the reference height [m]
- `scheme`: Surface flux solver scheme (default: `PointValueScheme()`)

# Formula:

    Cd = (κ / F_m)^2

where `F_m` is the dimensionless velocity profile.
"""
function drag_coefficient(
    param_set::APS,
    ζ,
    z0m,
    Δz,
    scheme = UF.PointValueScheme(),
)
    uf_params = SFP.uf_params(param_set)
    κ = SFP.von_karman_const(param_set)

    F_m = UF.dimensionless_profile(uf_params, Δz, ζ, z0m, UF.MomentumTransport(), scheme)
    Cd = (κ / F_m)^2
    return Cd
end

"""
    drag_coefficient(inputs::SurfaceFluxInputs, speed)

Compute the drag coefficient `Cd` from friction velocity (presumed to be in `inputs.ustar`) 
and effective wind speed (including any gustiness factors).
"""
function drag_coefficient(inputs::SurfaceFluxInputs, speed)
    ustar = inputs.ustar
    return (ustar / speed)^2
end

"""
    heat_exchange_coefficient(param_set, ζ, z0m, z0h, Δz, scheme)

Compute the heat exchange coefficient `Ch` for scalar exchange.

# Formula:

    Ch = κ^2 / (F_m * F_h),

where `F_m` and `F_h` are the dimensionless profiles for momentum and scalars. 
For the finite-volume case, this corresponds to the formulation in 
Nishizawa & Kitamura (2018), Eqs. 21 & 22 (with Pr_0 absorbed into F_h).

# Arguments
- `param_set`: Parameter set
- `ζ`: Stability parameter `ζ = Δz / L_MO`
- `z0m`: Roughness length for momentum [m]
- `z0h`: Roughness length for scalars (heat/moisture) [m]
- `Δz`: Height difference between the surface and the reference height [m]
- `scheme`: Surface flux solver scheme (default: `PointValueScheme()`)
"""
function heat_exchange_coefficient(
    param_set::APS,
    ζ,
    z0m,
    z0h,
    Δz,
    scheme = UF.PointValueScheme(),
)
    uf_params = SFP.uf_params(param_set)
    κ = SFP.von_karman_const(param_set)

    F_m = UF.dimensionless_profile(uf_params, Δz, ζ, z0m, UF.MomentumTransport(), scheme)
    F_h = UF.dimensionless_profile(uf_params, Δz, ζ, z0h, UF.HeatTransport(), scheme)

    Ch = κ^2 / (F_m * F_h)
    return Ch
end

"""
    heat_conductance(param_set, ζ, ustar, inputs, z0m, z0h, scheme)

Compute the heat conductance `g_h` (speed * Ch), including any gustiness factor in the wind speed.
Calculates windspeed and exchange coefficient internally from Monin-Obukhov variables.
"""
function heat_conductance(
    param_set::APS,
    ζ,
    ustar,
    inputs::SurfaceFluxInputs,
    z0m,
    z0h,
    scheme = UF.PointValueScheme(),
)
    # Compute Ch
    Ch = heat_exchange_coefficient(param_set, ζ, z0m, z0h, inputs.Δz, scheme)

    # Compute windspeed with gustiness (using windspeed helper which handles b_flux)
    current_speed = windspeed(param_set, ζ, ustar, inputs)

    return Ch * current_speed
end
