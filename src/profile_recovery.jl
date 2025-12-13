"""
    compute_profile_value(param_set, L_MO, z0, Δz, scale, val_sfc, transport)

Compute the value of a variable (momentum or scalar) at height `Δz` (height above surface).

# Arguments
- `param_set`: Parameter set
- `L_MO`: Monin-Obukhov length [m]
- `z0`: Roughness length [m]
- `Δz`: Height above the surface [m]
- `scale`: Similarity scale (u_star, theta_star, etc.)
- `val_sfc`: Surface value of the variable
- `transport`: Transport type (`MomentumTransport` or `HeatTransport`, the latter being 
    used for scalar transport)

# Formula:

    X(Δz) = (scale / κ) * F_z + val_sfc

where `F_z` is the dimensionless profile at height `Δz`.
"""
function compute_profile_value(
    param_set::APS,
    L_MO,
    z0,
    Δz,
    scale,
    val_sfc,
    transport,
    scheme = UF.PointValueScheme(),
)
    uf_params = SFP.uf_params(param_set)
    κ = SFP.von_karman_const(param_set)
    ζ = Δz / L_MO

    F = UF.dimensionless_profile(uf_params, Δz, ζ, z0, transport, scheme)

    return F * scale / κ + val_sfc
end
