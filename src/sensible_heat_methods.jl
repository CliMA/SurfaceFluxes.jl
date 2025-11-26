function sensible_heat_flux(
    param_set::APS,
    thermo_params,
    inputs::SurfaceFluxInputs,
    Ch,
    cp_m_in::FT,
    cp_m_sfc::FT,
    T_in::FT,
    T_sfc::FT,
    ρ_sfc::FT,
    gustiness::FT,
    E,
) where {FT}
    if inputs.shf !== nothing
        return inputs.shf
    end
    grav = SFP.grav(param_set)
    T_0 = SFP.T_0(param_set)
    LH_v0 = SFP.LH_v0(param_set)
    hv_sfc = TD.specific_enthalpy_vapor(thermo_params, T_sfc)
    ΔΦ = grav * inputs.Δz
    ΔDSE = cp_m_in * (T_in - T_0) - cp_m_sfc * (T_sfc - T_0) + ΔΦ
    Φ_sfc = surface_geopotential(inputs)
    return -ρ_sfc * Ch * windspeed(inputs, gustiness) * ΔDSE +
           (hv_sfc + Φ_sfc - LH_v0) * E
end
