"""
    compute_buoyancy_flux(param_set, shf, lhf, ts_in, ts_sfc, scheme)

Returns the buoyancy flux when the surface fluxes are known.
"""
function compute_buoyancy_flux(param_set, shf, lhf, ts_in, ts_sfc, scheme)
    thermo_params = SFP.thermodynamics_params(param_set)
    grav = SFP.grav(param_set)
    ε_vd = SFP.Rv_over_Rd(param_set)
    cp_m = TD.cp_m(thermo_params, ts_in)
    L_v = TD.latent_heat_vapor(thermo_params, ts_in)
    ρ_sfc = TD.air_density(thermo_params, ts_sfc)
    T_in = TD.air_temperature(thermo_params, ts_in)
    return grav / ρ_sfc * (shf / cp_m / T_in + (ε_vd - 1) * lhf / L_v)
end

function compute_buoyancy_flux(
    param_set,
    sc::Union{FluxesAndFrictionVelocity, Fluxes},
    scheme,
)
    return compute_buoyancy_flux(
        param_set,
        sc.shf,
        sc.lhf,
        ts_in(sc),
        ts_sfc(sc),
        scheme,
    )
end

