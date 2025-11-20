"""
    latent_heat_flux(param_set, Ch, sc, scheme)

In cases where surface fluxes are known,
return the known latent heat flux.
"""
function latent_heat_flux(
    param_set,
    Ch,
    sc::Union{Fluxes, FluxesAndFrictionVelocity},
    scheme,
)
    return sc.lhf
end

"""
    latent_heat_flux(param_set, Ch, sc, scheme)

Compute and return the latent heat flux
## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - Ch: Thermal exchange coefficient
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - scheme: Discretization scheme (currently supports FD and FV)
"""
function latent_heat_flux(
    param_set,
    Ch,
    sc::Union{ValuesOnly},
    scheme,
)
    LH_v0 = SFP.LH_v0(param_set)
    E = evaporation(param_set, sc, Ch)
    lhf = LH_v0 * E
    return lhf
end

function latent_heat_flux(
    param_set,
    Ch,
    sc::Union{Coefficients},
    scheme,
)
    LH_v0 = SFP.LH_v0(param_set)
    E = evaporation(param_set, sc, Ch)
    lhf = LH_v0 * E
    return lhf
end

