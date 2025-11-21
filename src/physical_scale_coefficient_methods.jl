"""
    compute_physical_scale_coeff(param_set, sc, L_MO, 𝓁, transport, ::LayerAverageScheme)

Computes the coefficient for the physical scale of a variable based on Nishizawa(2018)
for the FV scheme.

## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - L_MO: Monin-Obukhov length
  - 𝓁: Roughness length
  - transport: Transport type, (e.g. Momentum or Heat, used to determine physical scale coefficients)
  - scheme: Discretization scheme (currently supports FD and FV)
"""
function compute_physical_scale_coeff(
    param_set::APS,
    sc::Union{ValuesOnly, Fluxes, FluxesAndFrictionVelocity},
    L_MO,
    𝓁,
    transport,
    ::LayerAverageScheme,
)
    𝜅 = SFP.von_karman_const(param_set)
    uf = SFP.uf_params(param_set)
    π_group = UF.π_group(uf, transport)
    R_z0 = 1 - 𝓁 / Δz(sc)
    denom1 = log(Δz(sc) / 𝓁)
    denom2 = -UF.Psi(uf, Δz(sc) / L_MO, transport)
    denom3 =
        𝓁 / Δz(sc) *
        UF.Psi(uf, 𝓁 / L_MO, transport)
    denom4 = R_z0 * (UF.psi(uf, 𝓁 / L_MO, transport) - 1)
    Σterms = denom1 + denom2 + denom3 + denom4
    return 𝜅 / (π_group * Σterms)
end

"""
    compute_physical_scale_coeff(param_set, sc, L_MO, 𝓁, transport, ::PointValueScheme)

Computes the coefficient for the physical scale of a variable based on Byun (1990)
for the Finite Differences scheme.

## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - L_MO: Monin-Obukhov length
  - 𝓁: Roughness length
  - transport: Transport type, (e.g. Momentum or Heat, used to determine physical scale coefficients)
  - scheme: Discretization scheme (currently supports FD and FV)
"""
function compute_physical_scale_coeff(
    param_set,
    sc::Union{ValuesOnly, Fluxes, FluxesAndFrictionVelocity},
    L_MO,
    𝓁,
    transport,
    ::PointValueScheme,
)
    𝜅 = SFP.von_karman_const(param_set)
    FT = eltype(𝜅)
    uf = SFP.uf_params(param_set)
    π_group = UF.π_group(uf, transport)
    denom1 = log(FT(Δz(sc) / 𝓁))
    denom2 = -UF.psi(uf, FT(Δz(sc) / L_MO), transport)
    denom3 = UF.psi(uf, FT(𝓁 / L_MO), transport)
    Σterms = denom1 + denom2 + denom3
    return 𝜅 / (π_group * Σterms)
end
