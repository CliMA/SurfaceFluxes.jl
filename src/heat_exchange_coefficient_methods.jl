"""
    heat_exchange_coefficient(param_set, L_MO, u★, sc, scheme)

Compute and return Ch, the heat exchange coefficient.

Computes Ch from the friction velocity, heat scale, and wind speed using the Monin-Obukhov similarity theory.

## Arguments
- `param_set`: Abstract parameter set containing physical constants
- `L_MO`: Monin-Obukhov lengthscale
- `u★`: Friction velocity
- `sc`: Surface conditions container
- `scheme`: Discretization scheme (PointValueScheme or LayerAverageScheme)

## Returns
- `Ch`: Heat exchange coefficient
"""
function heat_exchange_coefficient(
    param_set,
    L_MO,
    u★,
    sc::Union{Fluxes, ValuesOnly, FluxesAndFrictionVelocity},
    scheme,
)
    thermo_params = SFP.thermodynamics_params(param_set)
    transport = UF.HeatTransport()
    𝜅 = SFP.von_karman_const(param_set)
    𝓁u = compute_z0(u★, param_set, sc, sc.roughness_model, UF.MomentumTransport())
    𝓁θ = compute_z0(u★, param_set, sc, sc.roughness_model, UF.HeatTransport())
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
    return Ch
end

"""
    heat_exchange_coefficient(param_set, L_MO, u★, sc::Coefficients, scheme)

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
)
    return sc.Ch
end
