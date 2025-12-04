@inline heat_conductance(inputs::SurfaceFluxInputs, Ch, gustiness) =
    Ch * windspeed(inputs, gustiness)

"""
    momentum_exchange_coefficient(param_set, L_MO, u★, 𝓁u, inputs, scheme, tol_neutral, gustiness, ΔDSEᵥ)

Compute and return `Cd`, the momentum exchange coefficient, for the current
similarity state. For neutral conditions (`abs(ΔDSEᵥ) <= tol_neutral`), uses the
logarithmic law of the wall; otherwise uses the diagnosed friction velocity.
"""
function momentum_exchange_coefficient(
    param_set::APS,
    L_MO,
    u★,
    𝓁u,
    inputs::SurfaceFluxInputs,
    scheme::SolverScheme,
    tol_neutral,
    gustiness::FT,
    ΔDSEᵥ_val::FT,
) where {FT}
    ΔU = windspeed(inputs, gustiness)
    Cd = (u★ / ΔU)^2
    return Cd
end

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
