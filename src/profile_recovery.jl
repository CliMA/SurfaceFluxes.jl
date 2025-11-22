
"""
    recover_profile(param_set, sc, L_MO, Z, X_in, X_sfc, transport, scheme)

Recover profiles of variable X given values of Z coordinates. Follows Nishizawa equation (21,22)
## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - L_MO: Monin-Obukhov length
  - Z: Z coordinate(s) (within surface layer) for which variable values are required
  - X_star: Scale parameter for variable X
  - X_sfc: For variable X, values at surface nodes
  - transport: Transport type, (e.g. Momentum or Heat, used to determine physical scale coefficients)
  - scheme: Discretization scheme (currently supports FD and FV)
"""
function recover_profile(
    param_set::APS,
    sc::AbstractSurfaceConditions,
    L_MO,
    𝓁,
    Z,
    X_star,
    X_sfc,
    transport,
    scheme::Union{LayerAverageScheme, PointValueScheme},
)
    uf = SFP.uf_params(param_set)
    𝜅 = SFP.von_karman_const(param_set)
    num1 = log(Z / 𝓁)
    num2 = -UF.psi(uf, Z / L_MO, transport)
    num3 = UF.psi(uf, 𝓁 / L_MO, transport)
    Σnum = num1 + num2 + num3
    return Σnum * X_star / 𝜅 + X_sfc
end
