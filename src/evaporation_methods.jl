"""
    evaporation(param_set, sc, Ch)

Compute and return the evaporation. When `sc` is `Fluxes` or `FluxesAndFrictionVelocity`,
evaporation is directly calculated from the latent heat flux.
## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - Ch: Thermal exchange coefficient
"""
function evaporation(
    param_set,
    sc::Union{Fluxes, FluxesAndFrictionVelocity},
    Ch,
)
    LH_v0 = SFP.LH_v0(param_set)
    return sc.lhf / LH_v0
end

"""
    evaporation(param_set, sc, Ch)

Compute and return the evaporation. When `sc` is `ValuesOnly` or `Coefficients`, a `beta` factor
is used to represent the resistance of the surface.
## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - Ch: Thermal exchange coefficient
"""
function evaporation(param_set, sc::Union{ValuesOnly}, Ch)
    thermo_params = SFP.thermodynamics_params(param_set)
    ρ_sfc = TD.air_density(thermo_params, ts_sfc(sc))
    return -ρ_sfc * Ch * windspeed(sc) * Δqt(param_set, sc) * sc.beta
end

function evaporation(param_set, sc::Union{Coefficients}, Ch)
    thermo_params = SFP.thermodynamics_params(param_set)
    ρ_sfc = TD.air_density(thermo_params, ts_sfc(sc))
    return -ρ_sfc * Ch * windspeed(sc) * Δqt(param_set, sc) * sc.beta
end
