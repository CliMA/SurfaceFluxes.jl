"""
    drag_coefficient(param_set, L_MO, z0m, Δz, scheme)

Compute the drag coefficient `Cd` for momentum exchange.

# Arguments
- `param_set`: Parameter set
- `L_MO`: Monin-Obukhov length [m]
- `z0m`: Roughness length for momentum [m]
- `Δz`: Height difference between the surface and the reference height [m]
- `scheme`: Surface flux solver scheme (default: `PointValueScheme()`)

# Formula:

    Cd = (κ / F_m)^2

where `F_m` is the dimensionless velocity profile.
"""
function drag_coefficient(
    param_set::APS,
    L_MO,
    z0m,
    Δz,
    scheme = UF.PointValueScheme(),
)
    uf_params = SFP.uf_params(param_set)
    κ = SFP.von_karman_const(param_set)
    ζ = Δz / L_MO

    F_m = UF.dimensionless_profile(uf_params, Δz, ζ, z0m, UF.MomentumTransport(), scheme)
    Cd = (κ / F_m)^2
    return Cd
end

"""
    heat_exchange_coefficient(param_set, L_MO, z0m, z0h, Δz, scheme)

Compute the heat exchange coefficient `Ch` for scalar exchange.

# Formula:

    Ch = κ^2 / (F_m * F_h),

where `F_m` and `F_h` are the dimensionless profiles for momentum and scalars. 
For the finite-volume case, this corresponds to the formulation in 
Nishizawa & Kitamura (2018), Eqs. 21 & 22 (with Pr_0 absorbed into F_h).

# Arguments
- `param_set`: Parameter set
- `L_MO`: Monin-Obukhov length [m]
- `z0m`: Roughness length for momentum [m]
- `z0h`: Roughness length for scalars (heat/moisture) [m]
- `Δz`: Height difference between the surface and the reference height [m]
- `scheme`: Surface flux solver scheme (default: `PointValueScheme()`)
"""
function heat_exchange_coefficient(
    param_set::APS,
    L_MO,
    z0m,
    z0h,
    Δz,
    scheme = UF.PointValueScheme(),
)
    uf_params = SFP.uf_params(param_set)
    κ = SFP.von_karman_const(param_set)
    ζ = Δz / L_MO

    F_m = UF.dimensionless_profile(uf_params, Δz, ζ, z0m, UF.MomentumTransport(), scheme)
    F_h = UF.dimensionless_profile(uf_params, Δz, ζ, z0h, UF.HeatTransport(), scheme)

    Ch = κ^2 / (F_m * F_h)
    return Ch
end

"""
    heat_conductance(inputs, Ch, gustiness)

Compute the heat conductance (speed * Ch).
"""
@inline heat_conductance(inputs::SurfaceFluxInputs, Ch, gustiness) =
    Ch * windspeed(inputs, gustiness)
