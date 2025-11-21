"""
    heat_exchange_coefficient(param_set, L_MO, sc, scheme, tol_neutral)

Compute and return Ch, the heat exchange coefficient given the
Monin-Obukhov lengthscale.
"""
function heat_exchange_coefficient(
    param_set,
    L_MO,
    u★,
    sc::Union{Fluxes, ValuesOnly, FluxesAndFrictionVelocity},
    scheme,
    tol_neutral,
)
    thermo_params = SFP.thermodynamics_params(param_set)
    transport = UF.HeatTransport()
    κ = SFP.von_karman_const(param_set)
    𝓁u = compute_z0(u★, param_set, sc, sc.roughness_model, UF.MomentumTransport())
    𝓁θ = compute_z0(u★, param_set, sc, sc.roughness_model, UF.HeatTransport())
    if abs(ΔDSEᵥ(param_set, sc)) <= tol_neutral
        Ch = κ^2 / (log(Δz(sc) / 𝓁θ) * log(Δz(sc) / 𝓁u))
    else
        ϕ_heat = compute_physical_scale_coeff(
            param_set,
            sc,
            L_MO,
            𝓁θ,
            transport,
            scheme,
        )
        ustar = compute_ustar(param_set, L_MO, 𝓁u, sc, scheme)
        Ch = ustar * ϕ_heat / windspeed(sc)
    end
    return Ch
end

"""
    heat_exchange_coefficient(param_set, L_MO, sc, scheme)

Return Ch, the heat exchange coefficient given the
Monin-Obukhov lengthscale.
"""
function heat_exchange_coefficient(
    param_set,
    L_MO,
    u★,
    sc::Coefficients,
    scheme,
    tol_neutral,
)
    return sc.Ch
end
