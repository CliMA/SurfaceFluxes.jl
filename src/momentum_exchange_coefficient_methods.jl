"""
    momentum_exchange_coefficient(param_set, L_MO, u★, sc, scheme, tol_neutral)

Compute and return Cd, the momentum exchange coefficient.

Computes Cd from the friction velocity and wind speed using the Monin-Obukhov similarity theory.

## Arguments
- `param_set`: Abstract parameter set containing physical constants
- `L_MO`: Monin-Obukhov lengthscale
- `u★`: Friction velocity
- `sc`: Surface conditions container
- `scheme`: Discretization scheme (PointValueScheme or LayerAverageScheme)
- `tol_neutral`: Tolerance for neutral stability detection (unused, kept for API compatibility)

## Returns
- `Cd`: Momentum exchange coefficient
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
    ustar = compute_ustar(param_set, L_MO, 𝓁, sc, scheme)
    Cd = ustar^2 / windspeed(sc)^2
    return Cd
end

"""
    momentum_exchange_coefficient(param_set, L_MO, u★, sc::Coefficients, scheme, tol_neutral)

Return the momentum exchange coefficient from the surface conditions.

When surface conditions are provided as exchange coefficients, this method
simply returns the pre-computed Cd value.
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
