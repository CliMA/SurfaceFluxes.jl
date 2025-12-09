"""
    compute_physical_scale_coeff(
        param_set::APS,
        inputs,
        L_MO,
        z0,
        transport,
        scheme::SolverScheme,
    )

Computes the coefficient for the physical scale of a variable.
Returns `ϕ` such that `scale = Δvalue * ϕ`.
For example, `u★ = ΔU * ϕ_m`.

This is computed as:
```math
ϕ = \\frac{κ}{F(z, L, z_0)}
```
where `F` is the dimensionless profile function from `UniversalFunctions`.
"""
function compute_physical_scale_coeff(
    param_set::APS,
    inputs,
    L_MO,
    z0,
    transport,
    scheme::SolverScheme,
)
    κ = SFP.von_karman_const(param_set)
    uf = SFP.uf_params(param_set)
    Δz_layer = Δz(inputs)
    ζ = Δz_layer / L_MO
    # We call dimensionless_profile from UniversalFunctions
    profile = UF.dimensionless_profile(
        uf,
        Δz_layer,
        ζ,
        z0,
        transport,
        scheme,
    )
    return κ / profile
end

"""
    compute_ustar(param_set, L_MO, z0, inputs, scheme, gustiness)

Return the friction velocity implied by the current Monin-Obukhov solution.
If a friction velocity is prescribed via `SurfaceFluxInputs`, it is returned
directly; otherwise it is recomputed from the similarity coefficients.
"""
function compute_ustar(
    param_set::APS,
    L_MO,
    z0,
    inputs::SurfaceFluxInputs,
    scheme::SolverScheme,
    gustiness,
)
    # Per-input ustar check
    if inputs.ustar !== nothing
        return inputs.ustar
    end
    
    ΔU = windspeed(inputs, gustiness)
    ϕ = compute_physical_scale_coeff(
        param_set,
        inputs,
        L_MO,
        z0,
        UF.MomentumTransport(),
        scheme,
    )
    return ΔU * ϕ
end
