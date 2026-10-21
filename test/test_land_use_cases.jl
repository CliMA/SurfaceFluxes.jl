# Test land use cases
#
# Tests specific coupled surface schemes:
# 1. Beta Model: Parameterizes moisture availability factor (β)
# 2. Canopy Model: Solves coupled energy/moisture balance with stomatal resistance
# 3. Solving for snow surface temperature that satisfies a flux balance
using SurfaceFluxes
using SurfaceFluxes.UniversalFunctions: BusingerParams
using SurfaceFluxes: Parameters
using ClimaParams
import Thermodynamics as TD
using Test

FT = Float64
param_set = Parameters.SurfaceFluxesParameters(FT, BusingerParams)
thermo_params = Parameters.thermodynamics_params(param_set)

# Mock inputs, which would be defined in the model calling SurfaceFluxes
T_int = FT(283.771);
q_tot_int = FT(0.004);
P_int = FT(92326.6);
ρ_int = TD.air_density(thermo_params, T_int, P_int);  # approximate, ignoring moisture
u_int = (FT(3), FT(1))
Δz = FT(18)
u_sfc = (FT(0), FT(0))
ρ_sfc = ρ_int;
Φ_sfc = FT(0)

@testset "Beta Model for Evaporation" begin
    d = FT(0)
    T_sfc_guess = T_int
    q_sat = TD.q_vap_saturation(
        thermo_params,
        T_sfc_guess,
        ρ_int,
        TD.Liquid(),
    )
    q_vap_sfc_guess = q_sat
    positional_default_args = (
        roughness_inputs = nothing,
        conf = SurfaceFluxes.default_surface_flux_config(eltype(param_set)),
        scheme = SurfaceFluxes.PointValueScheme(),
        solver_opts = nothing,
        flux_specs = nothing,
        update_T_sfc = nothing,
    )

    # Updates q_vap_sfc based on beta factor:
    # q_sfc = β * q_sat + (1 - β) * q_air
    function land_update_q_vap_sfc(ζ, param_set, thermo_params, inputs, β, q_sat)
        q_vap_int = inputs.q_tot_int - inputs.q_liq_int - inputs.q_ice_int
        q = β * q_sat + (1 - β) * q_vap_int
        return q
    end

    β = [FT(0), FT(0.5), FT(1)]
    evap = FT.(zeros(length(β)))
    for i in 1:length(β)
        update_q_vap_sfc(ζ, param_set, thermo_params, inputs, T_sfc, args...) =
            land_update_q_vap_sfc(ζ, param_set, thermo_params, inputs, β[i], q_sat)
        output = SurfaceFluxes.surface_fluxes(
            param_set,
            T_int,
            q_tot_int,
            FT(0),
            FT(0),
            ρ_int,
            T_sfc_guess,
            q_vap_sfc_guess,
            Φ_sfc,
            Δz,
            d,
            u_int,
            u_sfc,
            positional_default_args...,
            update_q_vap_sfc,
        )
        evap[i] = output.evaporation
    end

    # Compute evaporation directly with no update function (q = q_guess)
    output_potential = SurfaceFluxes.surface_fluxes(
        param_set,
        T_int,
        q_tot_int,
        FT(0),
        FT(0),
        ρ_int,
        T_sfc_guess,
        q_vap_sfc_guess,
        Φ_sfc,
        Δz,
        d,
        u_int,
        u_sfc,
    )
    @test evap[1] ≈ 0 # Δq = 0
    @test evap[2] < evap[3] # β  = 0.5 should reduce evaporation compared to β = 1
    @test evap[3] == output_potential.evaporation # using the beta function approach with β = 1 should agree with no beta function at all if q_sfc = q_sat
end

@testset "Canopy Model Evaporation" begin
    positional_default_args = (
        conf = SurfaceFluxes.default_surface_flux_config(eltype(param_set)),
        scheme = SurfaceFluxes.PointValueScheme(),
        solver_opts = nothing,
        flux_specs = nothing,
    )

    function land_update_T_sfc(
        ζ,
        param_set,
        thermo_params,
        inputs,
        scheme,
        u_star,
        z_0m,
        z_0h,
        leaf_Cd,
        LAI,
        T_canopy,
    )
        # Solve energy balance for T_sfc:
        # H + G = Rn (simplified here to coupling with canopy)
        # T_sfc is updated to balance heat fluxes given canopy conductance
        Φ_sfc = SurfaceFluxes.surface_geopotential(param_set, inputs)
        Φ_int = SurfaceFluxes.interior_geopotential(param_set, inputs)
        T_int = inputs.T_int
        g_h =
            SurfaceFluxes.heat_conductance(param_set, ζ, u_star, inputs, z_0m, z_0h, scheme)
        u = max(u_star, 1)
        g_land = leaf_Cd * u * LAI
        ΔΦ = Φ_int - Φ_sfc
        cp_d = TD.Parameters.cp_d(thermo_params)
        T_sfc = (T_int + T_canopy * g_land / g_h + ΔΦ / cp_d) / (1 + g_land / g_h)
        return T_sfc
    end

    function land_update_q_vap_sfc(
        ζ,
        param_set,
        thermo_params,
        inputs,
        scheme,
        T_sfc,
        u_star,
        z_0m,
        z_0h,
        leaf_Cd,
        LAI,
        r_stomata_canopy,
        q_canopy,
    )

        FT = eltype(param_set)
        u_star_safe = max(u_star, FT(1e-2))
        r_land = 1 / (leaf_Cd * u_star_safe) / max(LAI, eps(FT)) + r_stomata_canopy

        g_h =
            SurfaceFluxes.heat_conductance(param_set, ζ, u_star, inputs, z_0m, z_0h, scheme)

        # q_vap_int (atmosphere humidity)
        q_vap_int = inputs.q_tot_int

        # Solve for q_sfc analytically to satisfy balance of fluxes:
        # Flux_aero = ρ * g_h * (q_sfc - q_atm)
        # Flux_stom = ρ * (q_canopy - q_sfc) / r_land
        # Equating fluxes: g_h * (q_sfc - q_atm) = (q_canopy - q_sfc) / r_land
        # q_sfc * (g_h + 1/r_land) = q_canopy/r_land + g_h * q_atm
        # q_sfc = (q_canopy + g_h * r_land * q_atm) / (1 + g_h * r_land)

        q_new = (q_canopy + g_h * r_land * q_vap_int) / (1 + g_h * r_land)
        return q_new
    end

    r_stomata_canopy = FT(5)
    d = FT(0.67 * 10)
    z_0m = FT(0.13 * 10)
    z_0b = FT(0.1) * z_0m
    leaf_Cd = FT(0.03)
    roughness_inputs = SurfaceFluxes.ConstantRoughnessParams{FT}(z_0m, z_0b)
    T_canopy = T_int + FT(2)
    q_canopy = TD.q_vap_saturation(
        thermo_params,
        T_canopy,
        ρ_int,
        TD.Liquid(),
    )
    for LAI in [FT(2), FT(0)]
        update_T_sfc(ζ, param_set, thermo_params, inputs, scheme, u_star, z_0m, z_0h) =
            land_update_T_sfc(
                ζ,
                param_set,
                thermo_params,
                inputs,
                scheme,
                u_star,
                z_0m,
                z_0h,
                leaf_Cd,
                LAI,
                T_canopy,
            )
        update_q_vap_sfc(
            ζ,
            param_set,
            thermo_params,
            inputs,
            scheme,
            T_sfc,
            u_star,
            z_0m,
            z_0h,
        ) =
            land_update_q_vap_sfc(
                ζ,
                param_set,
                thermo_params,
                inputs,
                scheme,
                T_sfc,
                u_star,
                z_0m,
                z_0h,
                leaf_Cd,
                LAI,
                r_stomata_canopy,
                q_canopy,
            )

        output = SurfaceFluxes.surface_fluxes(
            param_set,
            T_int,
            q_tot_int,
            FT(0),
            FT(0),
            ρ_int,
            T_canopy,
            q_canopy,
            Φ_sfc,
            Δz,
            d,
            u_int,
            u_sfc,
            roughness_inputs,
            positional_default_args...,
            update_T_sfc,
            update_q_vap_sfc,
        )

        if LAI > 0
            @test output.shf > 0 # based on T_canopy, T_atmos
            @test output.lhf > 0 # based on q_canopy, q_atmos
        end

        if LAI ≈ FT(0) # Zero fluxes
            @test output.evaporation ≈ 0 atol = sqrt(eps(FT))
            @test abs.(output.shf) < sqrt(eps(FT))
        end
    end
end


@testset "Flux balance updates" begin
    function update_T_scheme(
        ζ,
        param_set,
        thermo_params,
        inputs,
        scheme,
        u_star,
        z_0m,
        z_0h,
        T̄,
        κ,
        d,
        σ,
        ϵ,
        SW_n,
        LW_d,
    )
        T_sfc = inputs.T_sfc_guess
        T_atmos = inputs.T_int
        ρ_atmos = inputs.ρ_int
        q_vap_atmos = inputs.q_tot_int - inputs.q_liq_int - inputs.q_ice_int
        q_sfc = TD.q_vap_saturation(thermo_params, T_sfc, ρ_atmos, FT(0), FT(0))
        ρ_sfc = SurfaceFluxes.surface_density(
            param_set,
            T_atmos,
            ρ_atmos,
            T_sfc,
            inputs.Δz,
            inputs.q_tot_int,
            inputs.q_liq_int,
            inputs.q_ice_int,
            q_sfc,
        )
        ∂q∂T = TD.∂q_vap_sat_∂T_from_L(
            thermo_params,
            q_sfc,
            TD.latent_heat_vapor(thermo_params, T_sfc),
            T_sfc,
        )

        g_h = SurfaceFluxes.heat_conductance(
            param_set, ζ, u_star, inputs, z_0m, z_0h, scheme,
        )
        E = SurfaceFluxes.evaporation(
            param_set,
            inputs,
            g_h,
            q_vap_atmos,
            q_sfc,
            ρ_sfc,
            inputs.moisture_model,
        )
        L = SurfaceFluxes.latent_heat_flux(
            param_set, inputs, E, inputs.moisture_model,
        )
        H = SurfaceFluxes.sensible_heat_flux(
            param_set, inputs, g_h, T_atmos, T_sfc, ρ_sfc, E,
        )

        _LH_v0 = TD.Parameters.LH_v0(thermo_params)
        cp_d = TD.Parameters.cp_d(thermo_params)
        ∂L∂T = ρ_sfc * g_h * _LH_v0 * ∂q∂T
        ∂H∂T = ρ_sfc * g_h * cp_d
        LW_n = ϵ * (LW_d - σ * T_sfc^4)
        ∂LW_n∂T = -4 * ϵ * σ * T_sfc^3

        # f(T) = L + H - R_n + κ (T_sfc - T̄) / d = 0
        ΔT =
            -(d * (-SW_n - LW_n + L + H) + κ * (T_sfc - T̄)) /
            (d * (-∂LW_n∂T + ∂L∂T + ∂H∂T) + κ)
        return T_sfc + ΔT
    end

    function update_q_scheme(
        ζ,
        param_set,
        thermo_params,
        inputs,
        scheme,
        T_sfc,
        u_star,
        z_0m,
        z_0h,
    )
        # Saturated surface humidity at the temperature returned by update_T_scheme
        return TD.q_vap_saturation(thermo_params, T_sfc, inputs.ρ_int, FT(0), FT(0))
    end

    positional_default_args = (
        conf = SurfaceFluxes.default_surface_flux_config(eltype(param_set)),
        scheme = SurfaceFluxes.PointValueScheme(),
        solver_opts = nothing,
        flux_specs = nothing,
    )
    d = FT(0.08)
    T̄ = T_int - FT(3)
    T_sfc_guess = T_int
    q_vap_sfc_guess = TD.q_vap_saturation(
        thermo_params,
        T_sfc_guess,
        ρ_int,
        FT(0),
        FT(0),
    )
    _σ = FT(5.67e-8) # Stefan–Boltzmann
    ϵ = FT(0.99)
    κ = FT(0.1)
    SW_n = FT(500)
    LW_d = FT(200)
    displ = FT(0.0)
    z_0m = FT(0.13)
    z_0b = FT(0.1) * z_0m
    roughness_inputs = SurfaceFluxes.ConstantRoughnessParams{FT}(z_0m, z_0b)
    update_T(args...) =
        update_T_scheme(args..., T̄, κ, d, _σ, ϵ, SW_n, LW_d)
    output = SurfaceFluxes.surface_fluxes(
        param_set,
        T_int,
        q_tot_int,
        FT(0),
        FT(0),
        ρ_int,
        T_sfc_guess,
        q_vap_sfc_guess,
        Φ_sfc,
        Δz,
        displ,
        u_int,
        u_sfc,
        roughness_inputs,
        positional_default_args...,
        update_T,
        update_q_scheme,
    )

    H = output.shf
    L = output.lhf
    LW_n = ϵ * (LW_d - _σ * output.T_sfc^4)
    F_int = -κ * (output.T_sfc - T̄) / d
    budget_residual = H + L - SW_n - LW_n - F_int
    budget_scale = max(abs(H + L), abs(SW_n + LW_n), abs(F_int), FT(1))
    @test abs(budget_residual) / budget_scale < FT(0.05)
end

@testset "Coupled soil skin energy and moisture balance under convective gustiness" begin
    # Low wind over warm/hot soil with dry-soil-layer vapor conductance g_liq:
    # the aerodynamic conductance varies strongly where the virtual potential
    # temperatures of the surface and the air cross.
    for FT_s in (Float32, Float64)
        ps = Parameters.SurfaceFluxesParameters(FT_s, BusingerParams)
        tp = Parameters.thermodynamics_params(ps)
        _σ = FT_s(5.670374419e-8)
        r = FT_s(0.016)
        ϵ = FT_s(0.96)
        LW_d = FT_s(450)
        P_atmos = FT_s(1e5)
        T_atmos = FT_s(308.27)
        q_atmos = FT_s(0.01956)
        ρ_atmos = TD.air_density(tp, T_atmos, P_atmos, q_atmos)
        Δz_s = FT_s(0.5)
        roughness = SurfaceFluxes.ConstantRoughnessParams{FT_s}(FT_s(0.01), FT_s(0.001))

        function soil_q_sat_and_deriv(T_s, inputs)
            ρ_s = SurfaceFluxes.surface_density(ps, inputs, T_s, inputs.q_vap_sfc_guess)
            q_s = TD.q_vap_saturation(tp, T_s, ρ_s, TD.Liquid())
            LH = TD.latent_heat_vapor(tp, T_s)
            ∂q∂T = TD.∂q_vap_sat_∂T_from_L(tp, q_s, LH, T_s)
            return ρ_s, q_s, ∂q∂T
        end

        function make_soil_callbacks(T_top, SW_n, g_liq)
            function update_T(ζ, p_set, t_params, inputs, scheme, u_star, z_0m, z_0b)
                T_s = inputs.T_sfc_guess
                q_air = SurfaceFluxes.interior_vapor_specific_humidity(inputs)
                ρ_s, q_s, ∂q∂T = soil_q_sat_and_deriv(T_s, inputs)
                g_h = SurfaceFluxes.heat_conductance(
                    p_set, ζ, u_star, inputs, z_0m, z_0b, scheme,
                )
                denom = g_liq + g_h
                β_e = ifelse(iszero(denom), one(denom), g_liq / denom)
                q_sfc = β_e * q_s + (1 - β_e) * q_air
                ∂q_sfc∂T = β_e * ∂q∂T
                E = SurfaceFluxes.evaporation(
                    p_set, inputs, g_h, q_air, q_sfc, ρ_s, inputs.moisture_model,
                )
                L = SurfaceFluxes.latent_heat_flux(p_set, inputs, E, inputs.moisture_model)
                H = SurfaceFluxes.sensible_heat_flux(
                    p_set, inputs, g_h, inputs.T_int, T_s, ρ_s, E,
                )
                LH_v0 = TD.Parameters.LH_v0(t_params)
                cp_m = TD.cp_m(t_params, q_sfc)
                ∂L∂T = ρ_s * g_h * LH_v0 * ∂q_sfc∂T
                ∂H∂T = ρ_s * g_h * cp_m + ∂q_sfc∂T * cp_m * (T_s - inputs.T_int)
                LW_n = ϵ * (LW_d - _σ * T_s^4)
                ∂LW_n∂T = -4 * ϵ * _σ * T_s^3
                f = SW_n - LW_n + L + H - (T_top - T_s) / r
                df = -∂LW_n∂T + ∂L∂T + ∂H∂T + 1 / r
                return T_s - f / df
            end
            function update_q(ζ, p_set, t_params, inputs, scheme, T_s, u_star, z_0m, z_0b)
                q_air = SurfaceFluxes.interior_vapor_specific_humidity(inputs)
                _, q_s, _ = soil_q_sat_and_deriv(T_s, inputs)
                g_h = SurfaceFluxes.heat_conductance(
                    p_set, ζ, u_star, inputs, z_0m, z_0b, scheme,
                )
                denom = g_liq + g_h
                β_e = ifelse(iszero(denom), one(denom), g_liq / denom)
                return β_e * q_s + (1 - β_e) * q_air
            end
            return update_T, update_q
        end

        for gust_min in (FT_s(0), FT_s(0.5)),
            T_top in (FT_s(300), FT_s(303), FT_s(304), FT_s(305), FT_s(308), FT_s(312.58)),
            u_spd in (FT_s(0), FT_s(0.02), FT_s(0.132), FT_s(0.5)),
            SW_n in (FT_s(-250), FT_s(0)),
            g_liq in (FT_s(1e-4), FT_s(0.01), FT_s(2.08))

            conf = SurfaceFluxes.SurfaceFluxConfig(
                roughness,
                SurfaceFluxes.FlooredDeardorffGustinessSpec(gust_min),
                SurfaceFluxes.MoistModel(),
                SurfaceFluxes.NoRoughnessSubLayer(),
                SurfaceFluxes.MaxHeatFluxStabilityCap(),
                SurfaceFluxes.ReferenceAboveSurface(),
            )
            update_T, update_q = make_soil_callbacks(T_top, SW_n, g_liq)
            q_guess = TD.q_vap_saturation(tp, T_top, ρ_atmos, TD.Liquid())
            out = SurfaceFluxes.surface_fluxes(
                ps,
                T_atmos, q_atmos, FT_s(0), FT_s(0), ρ_atmos,
                T_top, q_guess, FT_s(0), Δz_s, FT_s(0),
                (u_spd, FT_s(0)), (FT_s(0), FT_s(0)),
                nothing, conf, SurfaceFluxes.PointValueScheme(),
                nothing, nothing, update_T, update_q,
            )
            F = SW_n - ϵ * (LW_d - _σ * out.T_sfc^4) + out.lhf + out.shf
            G = (T_top - out.T_sfc) / r
            # Largest imbalances over this sweep: 1.0e-3 W/m² in Float64, and
            # 0.52 W/m² in Float32 at T_top = 303 K with zero wind and zero
            # gustiness floor, where Δθ_v at the solution is below the Float32
            # resolution of θ_v (≈ 3e-5 K).
            tol = FT_s === Float32 ? FT_s(1) : FT_s(2e-3)
            @test out.converged || out.ζ >= FT_s(100)
            @test isfinite(F) && abs(F - G) < tol
        end
    end
end
