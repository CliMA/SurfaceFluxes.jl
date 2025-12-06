"""
    compute_ustar(param_set, L_MO, 𝓁, inputs, scheme, gustiness)

Return the friction velocity implied by the current Monin-Obukhov solution.
If a friction velocity is prescribed via `SurfaceFluxInputs`, it is returned
directly; otherwise it is recomputed from the similarity coefficients.
"""
function compute_ustar(
    param_set::APS,
    L_MO,
    𝓁,
    inputs::SurfaceFluxInputs,
    scheme::SolverScheme,
    gustiness::FT,
) where {FT}
    if inputs.ustar !== nothing
        return inputs.ustar
    end
    ΔU = windspeed(inputs, gustiness)
    ϕ = compute_physical_scale_coeff(
        param_set,
        inputs,
        L_MO,
        𝓁,
        UF.MomentumTransport(),
        scheme,
    )
    return ΔU * ϕ
end
