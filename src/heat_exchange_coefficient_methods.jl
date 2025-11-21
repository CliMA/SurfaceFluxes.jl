"""
    heat_exchange_coefficient(param_set, L_MO, u★, sc, scheme, tol_neutral)

Compute and return Ch, the heat exchange coefficient.

For neutral conditions (when `abs(ΔDSEᵥ) <= tol_neutral`), uses the logarithmic
law of the wall. Otherwise, computes Ch from the friction velocity, heat scale,
and wind speed using the Monin-Obukhov similarity theory.

## Arguments
- `param_set`: Abstract parameter set containing physical constants
- `L_MO`: Monin-Obukhov lengthscale
- `u★`: Friction velocity
- `sc`: Surface conditions container
- `scheme`: Discretization scheme (PointValueScheme or LayerAverageScheme)
- `tol_neutral`: Tolerance for neutral stability detection (default: `cp_d / 100`)

## Returns
- `Ch`: Heat exchange coefficient
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
    𝜅 = SFP.von_karman_const(param_set)
    𝓁u = compute_z0(u★, param_set, sc, sc.roughness_model, UF.MomentumTransport())
    𝓁θ = compute_z0(u★, param_set, sc, sc.roughness_model, UF.HeatTransport())
    if abs(ΔDSEᵥ(param_set, sc)) <= tol_neutral
        Ch = 𝜅^2 / (log(Δz(sc) / 𝓁θ) * log(Δz(sc) / 𝓁u))
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
    heat_exchange_coefficient(param_set, L_MO, u★, sc::Coefficients, scheme, tol_neutral)

Return the heat exchange coefficient from the surface conditions.

When surface conditions are provided as exchange coefficients, this method
simply returns the pre-computed Ch value.
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
