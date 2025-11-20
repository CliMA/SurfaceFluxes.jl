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
