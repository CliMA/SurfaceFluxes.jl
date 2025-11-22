"""
    compute_ustar(
        param_set::AbstractSurfaceFluxesParameters,
        L_MO,
        𝓁,
        sc::AbstractSurfaceCondition,
        scheme,
    )

Return the friction velocity. This method is dispatched
by the surface condition:

## `sc::FluxesAndFrictionVelocity`

Friction velocity is known.

## `sc::Fluxes`

Compute given the Monin-Obukhov lengthscale.

## `sc::Coefficients`

Compute given the exchange coefficients.

## `sc::ValuesOnly`
Compute given the Monin-Obukhov lengthscale.
"""
function compute_ustar end

compute_ustar(param_set, L_MO, 𝓁, sc::FluxesAndFrictionVelocity, scheme) =
    sc.ustar

compute_ustar(param_set, L_MO, 𝓁, sc::Fluxes, scheme) =
    windspeed(sc) * compute_physical_scale_coeff(
        param_set,
        sc,
        L_MO,
        𝓁,
        UF.MomentumTransport(),
        scheme,
    )

compute_ustar(param_set, L_MO, 𝓁, sc::ValuesOnly, scheme) =
    windspeed(sc) * compute_physical_scale_coeff(
        param_set,
        sc,
        L_MO,
        𝓁,
        UF.MomentumTransport(),
        scheme,
    )

compute_ustar(param_set, L_MO, 𝓁, sc::Coefficients, scheme) =
    sqrt(sc.Cd) * (windspeed(sc))
