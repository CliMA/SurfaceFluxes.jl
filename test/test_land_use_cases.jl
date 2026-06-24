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
    q_vap_sfc_guess = TD.q_vap_saturation(
        thermo_params,
        T_sfc_guess,
        ρ_int,
        TD.Liquid(),
    )

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
    function land_update_q_vap_sfc(ζ, param_set, thermo_params, inputs, β)
        q_vap_int = inputs.q_tot_int - inputs.q_liq_int - inputs.q_ice_int
        q = β * inputs.q_vap_sfc_guess[1] + (1 - β) * q_vap_int # q_vap_sfc_guess is already the saturated value
        return q
    end

    β = [FT(0), FT(0.5), FT(1)]
    evap = FT.(zeros(length(β)))
    for i in 1:length(β)
        update_q_vap_sfc(ζ, param_set, thermo_params, inputs, T_sfc, args...) =
            land_update_q_vap_sfc(ζ, param_set, thermo_params, inputs, β[i])
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
        Φ_sfc = SurfaceFluxes.surface_geopotential(inputs)
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

    # What we did before for comparison:
    #    output_no_resistance = SurfaceFluxes.surface_fluxes(param_set, T_int, q_tot_int, ρ_int, T_canopy, q_canopy, Φ_sfc, Δz, d, u_int, u_sfc, roughness_inputs)
    #    r_e = r_stomata_canopy + 1/(leaf_Cd * max(output_no_resistance.ustar, 1))/max(LAI, eps(FT))
    #    r_h = 1/(leaf_Cd * max(output_no_resistance.ustar, 1))/max(LAI, eps(FT))
    #    r_ae = 1/output_no_resistance.Ch/max(1, sqrt(u_int[1]^2 + u_int[2]^2))
    #    pred_shf = output_no_resistance.shf *r_ae/(r_h+r_ae)
    #    pred_e = output_no_resistance.evaporation *r_ae/(r_e+r_ae)
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
        T_sfc = inputs.T_sfc_guess[1]
        T_atmos = inputs.T_int
        ρ_atmos = inputs.ρ_int
        q_atmos = inputs.q_tot_int
        Δz = inputs.Δz
        P_atmos =
            TD.air_pressure(thermo_params, T_atmos, ρ_atmos, q_atmos)
        q = TD.q_vap_saturation(
            thermo_params,
            T_sfc,
            ρ_atmos,
            FT(0),
            FT(0),
        )
        ∂q∂T = TD.∂q_vap_sat_∂T_from_L(
            thermo_params,
            q,
            TD.latent_heat_vapor(thermo_params, T_sfc),
            T_sfc,
        )

        g_h = SurfaceFluxes.heat_conductance(
            param_set,
            ζ,
            u_star,
            inputs,
            z_0m,
            z_0b,
            scheme,
        )
        b_flux = SurfaceFluxes.buoyancy_flux(param_set, ζ, u_star, inputs)
        E = SurfaceFluxes.evaporation(
            param_set,
            inputs,
            g_h,
            q_atmos,
            q_sfc,
            ρ_sfc,
            inputs.moisture_model,
        )
        L = SurfaceFluxes.latent_heat_flux(
            param_set,
            inputs,
            E,
            inputs.moisture_model,
        )
        H = SurfaceFluxes.sensible_heat_flux(
            param_set,
            inputs,
            g_h,
            T_atmos,
            T_sfc,
            ρ_sfc,
            E,
        )
        _LH_v0 = TD.Parameters.LH_v0(thermo_params)
        cp_d = TD.Parameters.cp_d(thermo_params)
        ∂L∂T = ρ_sfc * g_h * _LH_v0 * ∂q∂T
        ∂H∂T = ρ_sfc * g_h * cp_d
        LW_n = ϵ * (LW_d - σ * T_sfc^4)
        ∂LW_n∂T = -4 * ϵ * σ * T_sfc^3
        # f(T) = L + H - R_n + κ(T_sfc - T̄)/d = 0
        ΔT =
            -(d * (-SW_n - LW_n + L + H) + κ * (T_sfc - T̄)) /
            (d * (-∂LW_n∂T + ∂L∂T + ∂H∂T) + κ)
        @show T_sfc + ΔT
        @show ΔT
        inputs.T_sfc_guess .= [T_sfc + ΔT]
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
        q = TD.q_vap_saturation(
            thermo_params,
            inputs.T_sfc_guess[1], # this is the updated value after evaluated the update_T fn
            inputs.ρ_int,
            FT(0),
            FT(0),
        )
        return q
    end
    positional_default_args = (
        conf = SurfaceFluxes.default_surface_flux_config(eltype(param_set)),
        scheme = SurfaceFluxes.PointValueScheme(),
        solver_opts = nothing,
        flux_specs = nothing,
    )
    d = FT(0.08)
    T̄ = T_int - FT(3)
    T_sfc = T_int
    q_sfc = TD.q_vap_saturation(
        thermo_params,
        T_sfc,
        ρ_int,
        FT(0),
        FT(0),
    )
    _σ = FT(5.67e-8)# Stefan Boltzmann
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
        T_sfc,
        q_sfc,
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
    @show (H + L - LW_n - SW_n - F_int) / (H + L - LW_n - SW_n)
end
