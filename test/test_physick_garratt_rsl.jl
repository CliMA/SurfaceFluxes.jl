# Tests for the Physick & Garratt (1995) roughness sublayer (RSL) model.

module TestPhysickGarrattRSL

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions as UF
import ClimaParams as CP

@testset "Physick-Garratt RSL Model" begin
    FT = Float64
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)

    @testset "NoRoughnessSubLayer — zero correction" begin
        rsl = SF.NoRoughnessSubLayer()
        # Should return exactly zero for any inputs
        @test SF.rsl_profile_correction(rsl, FT(10), FT(0.1), UF.MomentumTransport()) == 0
        @test SF.rsl_profile_correction(rsl, FT(10), FT(0.1), UF.HeatTransport()) == 0
        @test SF.rsl_profile_correction(rsl, FT(100), FT(0.001), UF.MomentumTransport()) ==
              0
    end

    @testset "PhysickGarrattRSL — correction is negative (enhances exchange)" begin
        rsl = SF.PhysickGarrattRSL{FT}(c_m = FT(0.4), c_h = FT(0.4), z_RSL = FT(30.0))
        z0m = FT(0.5)   # momentum roughness length
        z0h = FT(0.05)  # scalar roughness length
        Δz_eff = FT(40.0)  # above RSL

        P_m = SF.rsl_profile_correction(rsl, Δz_eff, z0m, UF.MomentumTransport())
        P_h = SF.rsl_profile_correction(rsl, Δz_eff, z0h, UF.HeatTransport())

        # RSL correction must be negative (reduces F̂, enhances drag and scalar exchange)
        @test P_m < 0
        @test P_h < 0
    end

    @testset "PhysickGarrattRSL — correction saturates above RSL height" begin
        rsl = SF.PhysickGarrattRSL{FT}(c_m = FT(0.4), c_h = FT(0.4), z_RSL = FT(20.0))
        z0m = FT(0.5)

        P_at_rsl = SF.rsl_profile_correction(rsl, FT(20.0), z0m, UF.MomentumTransport())
        P_above1 = SF.rsl_profile_correction(rsl, FT(50.0), z0m, UF.MomentumTransport())
        P_above2 = SF.rsl_profile_correction(rsl, FT(100.0), z0m, UF.MomentumTransport())

        # Correction saturates at the RSL height
        @test P_at_rsl ≈ P_above1
        @test P_at_rsl ≈ P_above2
    end

    @testset "PhysickGarrattRSL — correction magnitude grows with c parameter" begin
        z0m = FT(0.5)
        Δz_eff = FT(40.0)
        z_RSL = FT(30.0)

        rsl_small = SF.PhysickGarrattRSL{FT}(c_m = FT(0.2), c_h = FT(0.2), z_RSL = z_RSL)
        rsl_large = SF.PhysickGarrattRSL{FT}(c_m = FT(0.6), c_h = FT(0.6), z_RSL = z_RSL)

        P_small = SF.rsl_profile_correction(rsl_small, Δz_eff, z0m, UF.MomentumTransport())
        P_large = SF.rsl_profile_correction(rsl_large, Δz_eff, z0m, UF.MomentumTransport())

        # Larger c → larger (more negative) correction
        @test P_large < P_small < 0
    end

    @testset "PhysickGarrattRSL — correction grows with z_RSL for fixed c" begin
        z0m = FT(0.5)
        Δz_eff = FT(100.0)  # well above both RSLs

        rsl_shallow =
            SF.PhysickGarrattRSL{FT}(c_m = FT(0.4), c_h = FT(0.4), z_RSL = FT(10.0))
        rsl_deep = SF.PhysickGarrattRSL{FT}(c_m = FT(0.4), c_h = FT(0.4), z_RSL = FT(40.0))

        P_shallow =
            SF.rsl_profile_correction(rsl_shallow, Δz_eff, z0m, UF.MomentumTransport())
        P_deep = SF.rsl_profile_correction(rsl_deep, Δz_eff, z0m, UF.MomentumTransport())

        # Deeper RSL → larger correction
        @test P_deep < P_shallow < 0
    end

    @testset "PhysickGarrattRSL — drag coefficient increases vs NoRSL" begin
        rsl = SF.PhysickGarrattRSL{FT}(c_m = FT(0.4), c_h = FT(0.4), z_RSL = FT(20.0))
        no_rsl = SF.NoRoughnessSubLayer()

        ζ = FT(0.0)   # neutral stability
        z0m = FT(0.5)
        Δz_eff = FT(30.0)
        scheme = UF.PointValueScheme()

        Cd_no_rsl = SF.drag_coefficient(param_set, ζ, z0m, Δz_eff, scheme, no_rsl)
        Cd_rsl = SF.drag_coefficient(param_set, ζ, z0m, Δz_eff, scheme, rsl)

        # RSL increases Cd
        @test Cd_rsl > Cd_no_rsl
    end

    @testset "PhysickGarrattRSL — heat exchange coefficient increases vs NoRSL" begin
        rsl = SF.PhysickGarrattRSL{FT}(c_m = FT(0.4), c_h = FT(0.4), z_RSL = FT(20.0))
        no_rsl = SF.NoRoughnessSubLayer()

        ζ = FT(0.0)
        z0m = FT(0.5)
        z0h = FT(0.05)
        Δz_eff = FT(30.0)
        scheme = UF.PointValueScheme()

        Ch_no_rsl =
            SF.heat_exchange_coefficient(param_set, ζ, z0m, z0h, Δz_eff, scheme, no_rsl)
        Ch_rsl = SF.heat_exchange_coefficient(param_set, ζ, z0m, z0h, Δz_eff, scheme, rsl)

        # RSL increases Ch
        @test Ch_rsl > Ch_no_rsl
    end

    @testset "PhysickGarrattRSL — full solve produces larger u* than NoRSL" begin
        # A forest-like canopy setup: h_c = 10m, z_RSL = 3 * h_c = 30m, d = 7m, z0m = 1m
        rsl_model = SF.PhysickGarrattRSL{FT}(c_m = FT(0.4), c_h = FT(0.4), z_RSL = FT(30.0))

        roughness = SF.ConstantRoughnessParams{FT}(z0m = FT(1.0), z0s = FT(0.1))
        gustiness = SF.ConstantGustinessSpec(FT(0.001))

        config_no_rsl = SF.SurfaceFluxConfig(roughness, gustiness, SF.DryModel())
        config_rsl = SF.SurfaceFluxConfig(roughness, gustiness, SF.DryModel(), rsl_model)

        # Typical forest-above-canopy measurement
        T_int = FT(290.0)
        T_sfc = FT(288.0)
        Δz = FT(40.0)   # measurement height above surface
        d = FT(7.0)    # displacement height
        u_int = (FT(5.0), FT(0.0))

        # Common inputs
        q = FT(0.0)
        ρ = FT(1.2)
        Φ_sfc = FT(0.0)

        result_no_rsl = SF.surface_fluxes(
            param_set, T_int, q, q, q, ρ, T_sfc, q, Φ_sfc, Δz, d, u_int, (FT(0), FT(0)),
            nothing, config_no_rsl,
        )
        result_rsl = SF.surface_fluxes(
            param_set, T_int, q, q, q, ρ, T_sfc, q, Φ_sfc, Δz, d, u_int, (FT(0), FT(0)),
            nothing, config_rsl,
        )

        # RSL should give larger u* (more drag)
        @test result_rsl.ustar > result_no_rsl.ustar
        # RSL should give larger Cd
        @test result_rsl.Cd > result_no_rsl.Cd
        # Both solves should converge
        @test result_rsl.converged
        @test result_no_rsl.converged
    end

    @testset "SurfaceFluxConfig backward-compatible constructors" begin
        roughness = SF.ConstantRoughnessParams{FT}()
        gustiness = SF.ConstantGustinessSpec(FT(1.0))

        # 2-arg constructor: should use MoistModel and NoRoughnessSubLayer defaults
        cfg2 = SF.SurfaceFluxConfig(roughness, gustiness)
        @test cfg2.moisture_model isa SF.MoistModel
        @test cfg2.rsl_model isa SF.NoRoughnessSubLayer

        # 3-arg constructor: explicit moisture model, NoRoughnessSubLayer default
        cfg3 = SF.SurfaceFluxConfig(roughness, gustiness, SF.DryModel())
        @test cfg3.moisture_model isa SF.DryModel
        @test cfg3.rsl_model isa SF.NoRoughnessSubLayer

        # 4-arg constructor: explicit RSL model
        rsl = SF.PhysickGarrattRSL{FT}(c_m = FT(0.4), c_h = FT(0.4), z_RSL = FT(20.0))
        cfg4 = SF.SurfaceFluxConfig(roughness, gustiness, SF.MoistModel(), rsl)
        @test cfg4.rsl_model isa SF.PhysickGarrattRSL
    end

    @testset "Float32 type stability" begin
        FT32 = Float32
        rsl =
            SF.PhysickGarrattRSL{FT32}(c_m = FT32(0.4), c_h = FT32(0.4), z_RSL = FT32(20.0))
        P = SF.rsl_profile_correction(rsl, FT32(30.0), FT32(0.5), UF.MomentumTransport())
        @test P isa FT32
        @test P < 0
    end
end

end # module
