"""
    heat_exchange_coefficient(param_set, L_MO, u★, 𝓁u, 𝓁θ, inputs, scheme, tol_neutral, gustiness, ΔDSEᵥ)

Compute and return `Ch`, the heat exchange coefficient, for the current
similarity state. Neutral and non-neutral regimes follow the log-law and MOST
formulations, respectively.
"""
function heat_exchange_coefficient(
    param_set::APS,
    L_MO,
    u★,
    𝓁u,
    𝓁θ,
    inputs::SurfaceFluxInputs,
    scheme::SolverScheme,
    tol_neutral,
    gustiness::FT,
    ΔDSEᵥ_val::FT,
) where {FT}
    transport = UF.HeatTransport()
    ΔU = windspeed(inputs, gustiness)
    ϕ_heat = compute_physical_scale_coeff(
        param_set,
        inputs,
        L_MO,
        𝓁θ,
        transport,
        scheme,
    )
    Ch = u★ * ϕ_heat / ΔU
    return Ch
end
