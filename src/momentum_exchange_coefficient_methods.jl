"""
    momentum_exchange_coefficient(param_set, L_MO, sc, scheme)

Compute and return Cd, the momentum exchange coefficient, given the
Monin-Obukhov lengthscale.
"""
function momentum_exchange_coefficient(
    param_set,
    L_MO,
    u★,
    sc::Union{Fluxes, ValuesOnly, FluxesAndFrictionVelocity},
    scheme,
    tol_neutral,
)
    thermo_params = SFP.thermodynamics_params(param_set)
    κ = SFP.von_karman_const(param_set)
    𝓁 = compute_z0(u★, param_set, sc, sc.roughness_model, UF.MomentumTransport())
    if abs(ΔDSEᵥ(param_set, sc)) <= tol_neutral
        Cd = (κ / log(Δz(sc) / 𝓁))^2
    else
        ustar = compute_ustar(param_set, L_MO, 𝓁, sc, scheme)
        Cd = ustar^2 / windspeed(sc)^2
    end
    return Cd
end

"""
    momentum_exchange_coefficient(param_set, L_MO, sc, scheme, tol_neutral)

Return Cd, the momentum exchange coefficient.
"""
function momentum_exchange_coefficient(
    param_set,
    L_MO,
    u★,
    sc::Coefficients,
    scheme,
    tol_neutral,
)
    return sc.Cd
end

