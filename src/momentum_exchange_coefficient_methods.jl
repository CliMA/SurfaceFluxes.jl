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
    κ = SFP.von_karman_const(param_set)
    ΔU = windspeed(inputs, gustiness)
    if abs(ΔDSEᵥ_val) <= tol_neutral
        Cd = (κ / log(inputs.Δz / 𝓁u))^2
    else
        Cd = (u★ / ΔU)^2
    end
    return Cd
end
