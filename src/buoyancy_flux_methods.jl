function compute_buoyancy_flux(
    param_set::APS,
    thermo_params,
    shf,
    lhf,
    cp_m_in,
    T_in,
    ρ_sfc,
)
    grav = SFP.grav(param_set)
    ε_vd = SFP.Rv_over_Rd(param_set)
    L_v = TD.latent_heat_vapor(thermo_params, T_in)
    return grav / ρ_sfc * (shf / (cp_m_in * T_in) + (ε_vd - 1) * lhf / L_v)
end
