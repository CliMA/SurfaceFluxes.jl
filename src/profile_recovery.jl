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
