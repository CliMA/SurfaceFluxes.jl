# Tests for the Harman & Finnigan (2007) roughness sublayer (RSL) model.

module TestHarmanFinniganRSL

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions as UF
import ClimaParams as CP

@testset "Harman-Finnigan RSL Model" begin
    FT = Float64
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)

    @testset "HarmanFinniganRSL — zero correction when c1=0" begin
        # When c1=0 the exponential is identity, so the integrand is zero everywhere.
        rsl = SF.HarmanFinniganRSL{FT}(c1_m = FT(0), c1_h = FT(0), z_RSL = FT(20.0))
        @test SF.rsl_profile_correction(rsl, FT(10), FT(0.1), UF.MomentumTransport()) == 0
        @test SF.rsl_profile_correction(rsl, FT(10), FT(0.1), UF.HeatTransport()) == 0
        @test SF.rsl_profile_correction(rsl, FT(100), FT(0.001), UF.MomentumTransport()) ==
              0
    end

    @testset "HarmanFinniganRSL — correction is negative (enhances exchange)" begin
        rsl = SF.HarmanFinniganRSL{FT}(c1_m = FT(0.5), c1_h = FT(0.5), z_RSL = FT(30.0))
        z0m = FT(0.5)
        z0h = FT(0.05)
        Δz_eff = FT(40.0)  # above RSL: correction has saturated

        P_m = SF.rsl_profile_correction(rsl, Δz_eff, z0m, UF.MomentumTransport())
        P_h = SF.rsl_profile_correction(rsl, Δz_eff, z0h, UF.HeatTransport())

        @test P_m < 0
        @test P_h < 0
    end

    @testset "HarmanFinniganRSL — correction saturates above RSL height" begin
        rsl = SF.HarmanFinniganRSL{FT}(c1_m = FT(0.5), c1_h = FT(0.5), z_RSL = FT(20.0))
        z0m = FT(0.5)

        P_at_rsl = SF.rsl_profile_correction(rsl, FT(20.0), z0m, UF.MomentumTransport())
        P_above1 = SF.rsl_profile_correction(rsl, FT(50.0), z0m, UF.MomentumTransport())
        P_above2 = SF.rsl_profile_correction(rsl, FT(100.0), z0m, UF.MomentumTransport())

        # Both use z_clip = z_RSL so the integral is identical
        @test P_at_rsl ≈ P_above1
        @test P_at_rsl ≈ P_above2
    end

    @testset "HarmanFinniganRSL — correction magnitude grows with c1" begin
        z0m = FT(0.5)
        Δz_eff = FT(40.0)
        z_RSL = FT(30.0)

        rsl_small = SF.HarmanFinniganRSL{FT}(c1_m = FT(0.2), c1_h = FT(0.2), z_RSL = z_RSL)
        rsl_large = SF.HarmanFinniganRSL{FT}(c1_m = FT(0.6), c1_h = FT(0.6), z_RSL = z_RSL)

        P_small = SF.rsl_profile_correction(rsl_small, Δz_eff, z0m, UF.MomentumTransport())
        P_large = SF.rsl_profile_correction(rsl_large, Δz_eff, z0m, UF.MomentumTransport())

        @test P_large < P_small < 0
    end

    @testset "HarmanFinniganRSL — weaker correction than PG for same parameter value" begin
        # Since exp(-x) > 1-x for x>0, the HF integrand [exp(-c(1-z/z*))-1]/z
        # is less negative than the PG integrand [-c(1-z/z*)]/z at every z.
        # Hence |P_HF| < |P_PG| for the same c1 = c parameter value.
        c = FT(0.5)
        z0m = FT(0.5)
        Δz_eff = FT(40.0)
        z_RSL = FT(30.0)

        rsl_hf = SF.HarmanFinniganRSL{FT}(c1_m = c, c1_h = c, z_RSL = z_RSL)
        rsl_pg = SF.PhysickGarrattRSL{FT}(c_m = c, c_h = c, z_RSL = z_RSL)

        P_hf = SF.rsl_profile_correction(rsl_hf, Δz_eff, z0m, UF.MomentumTransport())
        P_pg = SF.rsl_profile_correction(rsl_pg, Δz_eff, z0m, UF.MomentumTransport())

        @test P_pg < P_hf < 0  # both negative; PG more negative than HF
    end

    @testset "HarmanFinniganRSL — drag coefficient increases vs NoRSL" begin
        rsl = SF.HarmanFinniganRSL{FT}(c1_m = FT(0.5), c1_h = FT(0.5), z_RSL = FT(20.0))
        no_rsl = SF.NoRoughnessSubLayer()

        ζ = FT(0.0)
        z0m = FT(0.5)
        Δz_eff = FT(30.0)
        scheme = UF.PointValueScheme()

        Cd_no_rsl = SF.drag_coefficient(param_set, ζ, z0m, Δz_eff, scheme, no_rsl)
        Cd_rsl = SF.drag_coefficient(param_set, ζ, z0m, Δz_eff, scheme, rsl)

        @test Cd_rsl > Cd_no_rsl
    end

    @testset "HarmanFinniganRSL — heat exchange coefficient increases vs NoRSL" begin
        rsl = SF.HarmanFinniganRSL{FT}(c1_m = FT(0.5), c1_h = FT(0.5), z_RSL = FT(20.0))
        no_rsl = SF.NoRoughnessSubLayer()

        ζ = FT(0.0)
        z0m = FT(0.5)
        z0h = FT(0.05)
        Δz_eff = FT(30.0)
        scheme = UF.PointValueScheme()

        Ch_no_rsl =
            SF.heat_exchange_coefficient(param_set, ζ, z0m, z0h, Δz_eff, scheme, no_rsl)
        Ch_rsl = SF.heat_exchange_coefficient(param_set, ζ, z0m, z0h, Δz_eff, scheme, rsl)

        @test Ch_rsl > Ch_no_rsl
    end

    @testset "HarmanFinniganRSL — full solve produces larger u* and Cd than NoRSL" begin
        # Forest-like setup: h_c=10m, d=7m, z_RSL=10m above d
        rsl_model =
            SF.HarmanFinniganRSL{FT}(c1_m = FT(0.5), c1_h = FT(0.5), z_RSL = FT(10.0))

        roughness = SF.ConstantRoughnessParams{FT}(z0m = FT(1.0), z0s = FT(0.1))
        gustiness = SF.ConstantGustinessSpec(FT(0.001))

        config_no_rsl = SF.SurfaceFluxConfig(roughness, gustiness, SF.DryModel())
        config_rsl = SF.SurfaceFluxConfig(roughness, gustiness, SF.DryModel(), rsl_model)

        T_int = FT(290.0)
        T_sfc = FT(288.0)
        Δz = FT(40.0)
        d = FT(7.0)
        u_int = (FT(5.0), FT(0.0))
        q = FT(0.0);
        ρ = FT(1.2);
        Φ_sfc = FT(0.0)

        result_no_rsl = SF.surface_fluxes(
            param_set, T_int, q, q, q, ρ, T_sfc, q, Φ_sfc, Δz, d, u_int, (FT(0), FT(0)),
            nothing, config_no_rsl,
        )
        result_rsl = SF.surface_fluxes(
            param_set, T_int, q, q, q, ρ, T_sfc, q, Φ_sfc, Δz, d, u_int, (FT(0), FT(0)),
            nothing, config_rsl,
        )

        @test result_rsl.ustar > result_no_rsl.ustar
        @test result_rsl.Cd > result_no_rsl.Cd
        @test result_rsl.converged
        @test result_no_rsl.converged
    end

    @testset "HarmanFinniganRSL — SurfaceFluxConfig constructor" begin
        roughness = SF.ConstantRoughnessParams{FT}()
        gustiness = SF.ConstantGustinessSpec(FT(1.0))
        rsl = SF.HarmanFinniganRSL{FT}(c1_m = FT(0.5), c1_h = FT(0.5), z_RSL = FT(15.0))
        cfg = SF.SurfaceFluxConfig(roughness, gustiness, SF.DryModel(), rsl)
        @test cfg.rsl_model isa SF.HarmanFinniganRSL
    end

    @testset "HarmanFinniganRSL — Float32 type stability" begin
        FT32 = Float32
        rsl = SF.HarmanFinniganRSL{FT32}(
            c1_m = FT32(0.5),
            c1_h = FT32(0.5),
            z_RSL = FT32(20.0),
        )
        P = SF.rsl_profile_correction(rsl, FT32(30.0), FT32(0.5), UF.MomentumTransport())
        @test P isa FT32
        @test P < 0
    end
end

end # module
