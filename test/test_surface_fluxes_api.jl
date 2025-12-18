module TestSurfaceFluxesAPI

using Test
using SurfaceFluxes
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD
import Thermodynamics.Parameters as TDP
using ClimaParams

# Shared Parameter Set
FT = Float64
const param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)

@testset "Surface Flux Primitives API" begin

    T_int = FT(300)
    q_tot_int = FT(0.01)
    T_sfc_guess = FT(298)
    q_vap_sfc_guess = FT(0.012)
    ρ_int = FT(1.2) # Simplified density

    base_result = SF.surface_fluxes(
        param_set,
        T_int,
        q_tot_int,
        ρ_int,
        T_sfc_guess,
        q_vap_sfc_guess,
        FT(0),
        FT(10),
        FT(0),
    )

    @test base_result isa SF.SurfaceFluxConditions{FT}
    @test base_result.Cd > zero(FT)

    Ts_calls = Ref(0)
    qs_calls = Ref(0)
    update_T_sfc = function (ζ, param_set, thermo_params, inputs, scheme, u_star, z0m, z0s)
        Ts_calls[] += 1
        return inputs.T_int - FT(5)
    end
    update_q_vap_sfc =
        function (ζ, param_set, thermo_params, inputs, scheme, T_sfc, u_star, z0m, z0s)
            qs_calls[] += 1
            return inputs.q_tot_int + FT(2e-3)
        end

    hooked_result = SF.surface_fluxes(
        param_set,
        T_int,
        q_tot_int,
        ρ_int,
        T_sfc_guess,
        q_vap_sfc_guess,
        FT(0),
        FT(10),
        FT(0),
        (0, 0),
        (0, 0),
        nothing,
        SF.default_surface_flux_config(FT),
        SF.PointValueScheme(),
        nothing,
        nothing,
        update_T_sfc,
        update_q_vap_sfc,
    )

    @test Ts_calls[] > 0
    @test qs_calls[] > 0
    @test hooked_result.shf != base_result.shf

    config = SF.SurfaceFluxConfig(
        SF.COARE3RoughnessParams{FT}(),
        SF.gustiness_constant(FT(0.5)),
    )

    config_result = SF.surface_fluxes(
        param_set,
        T_int,
        q_tot_int,
        ρ_int,
        T_sfc_guess,
        q_vap_sfc_guess,
        FT(0),
        FT(10),
        FT(0),
        (0, 0),
        (0, 0),
        nothing,
        config,
    )

    @test config_result.Cd != base_result.Cd
end

@testset "DryModel" begin
    # Setup inputs
    T_int = 300.0
    q_tot_int = 0.0 # Dry atmosphere
    ρ_int = 1.2
    T_sfc_guess = 305.0
    q_vap_sfc_guess = 0.0
    Φ_sfc = 0.0
    Δz = 10.0
    d = 0.0
    u_int = (10.0, 0.0)
    u_sfc = (0.0, 0.0)

    # Configure with DryModel
    config = SF.SurfaceFluxConfig(
        SF.ConstantRoughnessParams(1e-3, 1e-3),
        SF.ConstantGustinessSpec(1.0),
        SF.DryModel(),
    )

    # Call surface_fluxes
    sf = SF.surface_fluxes(
        param_set,
        T_int, q_tot_int, ρ_int, T_sfc_guess, q_vap_sfc_guess, Φ_sfc, Δz, d,
        u_int, u_sfc, nothing, config,
    )

    @test sf.lhf == 0.0
    @test sf.evaporation == 0.0
    @test sf.shf != 0.0 # Should be non-zero as T diff exists
end

@testset "Prescribed Fluxes" begin
    # Fully prescribed
    shf_pre = 100.0
    lhf_pre = 200.0
    ustar_pre = 0.5

    flux_specs = SF.FluxSpecs{FT}(; shf = shf_pre, lhf = lhf_pre, ustar = ustar_pre)

    T_int = 300.0
    q_tot_int = 0.01
    ρ_int = 1.2
    T_sfc_guess = 300.0
    q_vap_sfc_guess = 0.01
    Φ_sfc = 0.0
    Δz = 10.0
    d = 0.0
    u_int = (10.0, 0.0)
    u_sfc = (0.0, 0.0)

    config = SF.default_surface_flux_config(FT)

    sf = SF.surface_fluxes(
        param_set,
        T_int, q_tot_int, ρ_int, T_sfc_guess, q_vap_sfc_guess, Φ_sfc, Δz, d,
        u_int, u_sfc, nothing, config,
        SF.PointValueScheme(), nothing, flux_specs,
    )

    @test sf.shf == shf_pre
    @test sf.lhf == lhf_pre
    @test sf.ustar == ustar_pre
    @test sf.Cd > 0
    @test sf.L_MO != 0
end

@testset "Prescribed Fluxes and Coefficients Mixed" begin
    # Prescribed Coefficients ONLY (Legacy path)
    Cd_pre = 1e-3
    Ch_pre = 1e-3
    flux_specs = SF.FluxSpecs{FT}(; Cd = Cd_pre, Ch = Ch_pre)

    T_int = 300.0
    q_tot_int = 0.01
    ρ_int = 1.2
    T_sfc_guess = 305.0
    q_vap_sfc_guess = 0.02
    Φ_sfc = 0.0
    Δz = 10.0
    d = 0.0
    u_int = (10.0, 0.0)
    u_sfc = (0.0, 0.0)

    config = SF.default_surface_flux_config(FT)

    sf = SF.surface_fluxes(
        param_set,
        T_int, q_tot_int, ρ_int, T_sfc_guess, q_vap_sfc_guess, Φ_sfc, Δz, d,
        u_int, u_sfc, nothing, config,
        SF.PointValueScheme(), nothing, flux_specs,
    )

    @test sf.Cd == Cd_pre
    @test sf.Ch == Ch_pre
    @test sf.shf != 0.0 # calculated
end

@testset "Prescribed Ustar Only" begin
    ustar_pre = 0.5
    flux_specs = SF.FluxSpecs{FT}(; ustar = ustar_pre)

    T_int = 300.0
    q_tot_int = 0.01
    ρ_int = 1.2
    T_sfc_guess = 305.0
    q_vap_sfc_guess = 0.02
    Φ_sfc = 0.0
    Δz = 10.0
    d = 0.0
    u_int = (10.0, 0.0)
    u_sfc = (0.0, 0.0)

    config = SF.default_surface_flux_config(FT)

    sf = SF.surface_fluxes(
        param_set,
        T_int, q_tot_int, ρ_int, T_sfc_guess, q_vap_sfc_guess, Φ_sfc, Δz, d,
        u_int, u_sfc, nothing, config,
        SF.PointValueScheme(), nothing, flux_specs,
    )

    @test sf.ustar == ustar_pre
    @test sf.Cd ≈ (ustar_pre / sqrt(10.0^2 + 1.0^2)) atol = 0.1 # Approx check, gustiness is 1.0
    @test sf.shf != 0

    # Test sensible_heat_flux default E
    # Manually call sensible_heat_flux to verify default arg
    shf_default_E = SF.sensible_heat_flux(
        param_set,
        SF.build_surface_flux_inputs(
            T_int, q_tot_int, zero(FT), zero(FT), ρ_int, T_sfc_guess, q_vap_sfc_guess,
            Φ_sfc, Δz, d,
            u_int, u_sfc, config,
            nothing, SF.FluxSpecs{FT}(), nothing, nothing,
        ),
        FT(0.01), # g_h
        T_int,
        T_sfc_guess,
        ρ_int,
        # E defaulted
    )
    @test shf_default_E isa FT

    # Test helper methods
    ζ_test = FT(0.1)
    ustar_test = FT(0.3)
    z0m_test = FT(1e-3)
    z0h_test = FT(1e-3)

    # Rebuild inputs for the test call since sf is SurfaceFluxConditions
    inputs_test = SF.build_surface_flux_inputs(
        T_int, q_tot_int, zero(FT), zero(FT), ρ_int, T_sfc_guess, q_vap_sfc_guess,
        Φ_sfc, Δz, d, u_int,
        u_sfc, config,
        nothing, SF.FluxSpecs{FT}(), nothing, nothing,
    )

    shf_helper = SF.sensible_heat_flux(
        param_set,
        ζ_test,
        ustar_test,
        inputs_test,
        z0m_test,
        z0h_test,
        T_sfc_guess,
        q_vap_sfc_guess,
        ρ_int,
        SF.PointValueScheme(),
    )
    @test shf_helper isa FT
    @test shf_helper isa FT
    @test shf_helper != 0

    # Test windspeed helper
    ws_helper = SF.windspeed(
        param_set,
        ζ_test,
        ustar_test,
        inputs_test,
    )
    @test ws_helper isa FT
    @test ws_helper > 0

    # Test heat_conductance helper
    gh_helper = SF.heat_conductance(
        param_set,
        ζ_test,
        ustar_test,
        inputs_test,
        z0m_test,
        z0h_test,
        SF.PointValueScheme(),
    )
    @test gh_helper isa FT
    @test gh_helper isa FT
    @test gh_helper > 0

    # Test evaporation helper
    evap_helper = SF.evaporation(
        param_set,
        ζ_test,
        ustar_test,
        inputs_test,
        z0m_test,
        z0h_test,
        q_vap_sfc_guess,
        ρ_int,
        SF.PointValueScheme(),
    )
    @test evap_helper isa FT
    # Evaporation might be positive or negative depending on gradient, just check it runs
    @test evap_helper isa FT
    # Evaporation might be positive or negative depending on gradient, just check it runs
    @test !isnan(evap_helper)

    # Test latent_heat_flux helper
    lhf_helper = SF.latent_heat_flux(
        param_set,
        ζ_test,
        ustar_test,
        inputs_test,
        z0m_test,
        z0h_test,
        q_vap_sfc_guess,
        ρ_int,
        SF.PointValueScheme(),
    )
    @test lhf_helper isa FT
    @test !isnan(lhf_helper)
end

end # module
