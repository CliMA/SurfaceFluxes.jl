function obukhov_similarity_solution(
    param_set,
    sc::Coefficients,
    scheme,
    args...,
)
    lhf = latent_heat_flux(param_set, sc.Ch, sc, scheme)
    shf = sensible_heat_flux(param_set, sc.Ch, sc, scheme)
    ustar = sqrt(sc.Cd) * windspeed(sc)
    buoyancy_flux = compute_buoyancy_flux(
        param_set,
        shf,
        lhf,
        ts_in(sc),
        ts_sfc(sc),
        scheme,
    )
    return (L★ = -ustar^3 / SFP.von_karman_const(param_set) / non_zero(buoyancy_flux), u★ = ustar)
end

compute_ustar(param_set, L_MO, 𝓁, sc::Coefficients, scheme) =
    sqrt(sc.Cd) * (windspeed(sc))
"""
    momentum_exchange_coefficient(param_set, L_MO, sc, scheme, tol_neutral)

Return Cd, the momentum exchange coefficient.
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
"""
    heat_exchange_coefficient(param_set, L_MO, sc, scheme)

Return Ch, the heat exchange coefficient given the
Monin-Obukhov lengthscale.
"""
function heat_exchange_coefficient(
    param_set,
    L_MO,
    u★,
    sc::Coefficients,
    scheme,
    tol_neutral,
)
    return sc.Ch
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
function evaporation(param_set, sc::Union{Coefficients}, Ch)
    thermo_params = SFP.thermodynamics_params(param_set)
    ρ_sfc = TD.air_density(thermo_params, ts_sfc(sc))
    return -ρ_sfc * Ch * windspeed(sc) * Δqt(param_set, sc) * sc.beta
end
function surface_conditions(
    param_set::APS{FT},
    sc::Coefficients,
    scheme::SolverScheme = PointValueScheme();
    tol_neutral = SFP.cp_d(param_set) / 100,
    tol = sqrt(eps(FT)),
    maxiter::Int = 10,
) where {FT}
    uft = SFP.uf_params(param_set)
    ustar = sqrt(sc.Cd) * windspeed(sc)
    X★ = obukhov_similarity_solution(param_set, sc, uft, scheme)
    Cd = momentum_exchange_coefficient(
        param_set,
        nothing,
        nothing,
        sc,
        scheme,
        tol_neutral,
    )
    Ch =
        heat_exchange_coefficient(
            param_set,
            nothing,
            nothing,
            sc,
            scheme,
            tol_neutral,
        )
    shf = sensible_heat_flux(param_set, Ch, sc, scheme)
    lhf = latent_heat_flux(param_set, Ch, sc, scheme)
    buoy_flux = compute_buoyancy_flux(
        param_set,
        shf,
        lhf,
        ts_in(sc),
        ts_sfc(sc),
        scheme,
    )
    ρτxz, ρτyz = momentum_fluxes(param_set, Cd, sc, scheme)
    E = evaporation(param_set, sc, Ch)
    return SurfaceFluxConditions(
        X★.L★,
        shf,
        lhf,
        buoy_flux,
        ρτxz,
        ρτyz,
        X★.u★,
        Cd,
        Ch,
        E,
    )
end
