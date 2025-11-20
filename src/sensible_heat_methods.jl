"""
    sensible_heat_flux(param_set, Ch, sc, scheme)

Compute and return the sensible heat fluxes
## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - Ch: Thermal exchange coefficient
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - scheme: Discretization scheme (currently supports FD and FV)
"""
function sensible_heat_flux(
    param_set,
    Ch,
    sc::Union{ValuesOnly},
    scheme,
)
    thermo_params = SFP.thermodynamics_params(param_set)
    grav = SFP.grav(param_set)
    T_0 = SFP.T_0(param_set)
    LH_v0 = SFP.LH_v0(param_set)
    cp_m_in = TD.cp_m(thermo_params, ts_in(sc))
    cp_m_sfc = TD.cp_m(thermo_params, ts_sfc(sc))
    T_in = TD.air_temperature(thermo_params, ts_in(sc))
    T_sfc = TD.air_temperature(thermo_params, ts_sfc(sc))
    ρ_sfc = TD.air_density(thermo_params, ts_sfc(sc))
    hv_sfc = TD.specific_enthalpy_vapor(thermo_params, T_sfc)
    ΔΦ = grav * Δz(sc)
    ΔDSE = cp_m_in * (T_in - T_0) - cp_m_sfc * (T_sfc - T_0) + ΔΦ
    Φ_sfc = grav * z_sfc(sc)
    E = evaporation(param_set, sc, Ch)
    return -ρ_sfc * Ch * windspeed(sc) * ΔDSE + (hv_sfc + Φ_sfc - LH_v0) * E
end

"""
    sensible_heat_flux(param_set, Ch, sc, scheme)

In cases where surface fluxes are known,
return the known sensible heat flux.
"""
function sensible_heat_flux(
    param_set,
    Ch,
    sc::Union{Fluxes, FluxesAndFrictionVelocity},
    scheme,
)
    return sc.shf
end

function sensible_heat_flux(
    param_set,
    Ch,
    sc::Union{Coefficients},
    scheme,
)
    thermo_params = SFP.thermodynamics_params(param_set)
    grav = SFP.grav(param_set)
    T_0 = SFP.T_0(param_set)
    LH_v0 = SFP.LH_v0(param_set)
    cp_m_in = TD.cp_m(thermo_params, ts_in(sc))
    cp_m_sfc = TD.cp_m(thermo_params, ts_sfc(sc))
    ρ_sfc = TD.air_density(thermo_params, ts_sfc(sc))
    T_in = TD.air_temperature(thermo_params, ts_in(sc))
    T_sfc = TD.air_temperature(thermo_params, ts_sfc(sc))
    hv_sfc = TD.specific_enthalpy_vapor(thermo_params, T_sfc)
    ΔΦ = grav * Δz(sc)
    ΔDSE = cp_m_in * (T_in - T_0) - cp_m_sfc * (T_sfc - T_0) + ΔΦ
    Φ_sfc = grav * z_sfc(sc)
    E = evaporation(param_set, sc, Ch)
    return -ρ_sfc * Ch * windspeed(sc) * ΔDSE + (hv_sfc + Φ_sfc - LH_v0) * E
end

