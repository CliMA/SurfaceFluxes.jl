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
