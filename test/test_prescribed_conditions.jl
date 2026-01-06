using Test

import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import ClimaParams as CP

FT = Float32
param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)

@testset "Prescribed Fluxes and Coefficients" begin
    # Create a base case (standard conditions)
    # We use this to compute reference "truth" values for ustar, fluxes, coefficients
    # and then feed them back in to verify the prescribed mode logic.

    T_sfc = FT(300.0)
    T_int = FT(300.0) - FT(1.0) # Unstable
    q_sfc = FT(0.015)
    q_int = FT(0.012)
    u_int = (FT(10.0), FT(0.0))
    u_sfc = (FT(0.0), FT(0.0))
    Δz = FT(10.0)
    ρ_int = FT(1.2)

    config = SF.SurfaceFluxConfig(
        SF.ConstantRoughnessParams(FT(1e-3), FT(1e-3)),
        SF.ConstantGustinessSpec(FT(1.0)),
        SF.MoistModel(),
    )

    options = SF.SolverOptions{FT}(maxiter = 100, tol = FT(1e-8))

    # 1. Run Standard MOST helper (Ground Truth)
    base_result = SF.surface_fluxes(
        param_set,
        T_int, q_int, FT(0), FT(0), ρ_int,
        T_sfc, q_sfc,
        FT(0), Δz, FT(0),
        u_int, u_sfc,
        nothing,
        config,
        SF.PointValueScheme(),
        options,
    )

    @test base_result.converged
    @test base_result.T_sfc == T_sfc
    @test base_result.q_vap_sfc == q_sfc

    # Note: "Flux-prescribed" (given shf/lhf, solve ustar) is not supported in current API.
    # We verify only "Flux+Ustar prescribed" and "Coefficient prescribed".

    @testset "Flux + Ustar Prescribed" begin
        # Pass shf, lhf, ustar via FluxSpecs
        specs = SF.FluxSpecs{FT}(
            shf = base_result.shf,
            lhf = base_result.lhf,
            ustar = base_result.ustar,
        )

        result = SF.surface_fluxes(
            param_set,
            T_int, q_int, FT(0), FT(0), ρ_int,
            T_sfc, q_sfc,
            FT(0), Δz, FT(0),
            u_int, u_sfc,
            nothing,
            config,
            SF.PointValueScheme(),
            nothing, # No solver opts needed
            specs,
        )

        # Verify strictly that inputs are passed through
        @test result.shf == base_result.shf
        @test result.lhf == base_result.lhf
        @test result.ustar == base_result.ustar

    end

    @testset "Coefficients Prescribed" begin
        # Pass Cd, Ch via FluxSpecs
        specs = SF.FluxSpecs{FT}(
            Cd = base_result.Cd,
            Ch = base_result.Ch,
        )

        result = SF.surface_fluxes(
            param_set,
            T_int, q_int, FT(0), FT(0), ρ_int,
            T_sfc, q_sfc,
            FT(0), Δz, FT(0),
            u_int, u_sfc,
            nothing,
            config,
            SF.PointValueScheme(),
            nothing,
            specs,
        )

        @test result.Cd == base_result.Cd
        @test result.Ch == base_result.Ch
        # Check derived fluxes
        @test isapprox(result.shf, base_result.shf; rtol = 1e-4)
        @test isapprox(result.lhf, base_result.lhf; rtol = 1e-4)
        @test isapprox(result.ustar, base_result.ustar; rtol = 1e-4)
    end

    @testset "Prescribed SHF, LHF, and Cd" begin
        # Pass shf, lhf, Cd via FluxSpecs
        specs = SF.FluxSpecs{FT}(
            shf = base_result.shf,
            lhf = base_result.lhf,
            Cd = base_result.Cd,
        )

        result = SF.surface_fluxes(
            param_set,
            T_int, q_int, FT(0), FT(0), ρ_int,
            T_sfc, q_sfc,
            FT(0), Δz, FT(0),
            u_int, u_sfc,
            nothing,
            config,
            SF.PointValueScheme(),
            nothing,
            specs,
        )

        # Inputs pass through
        @test result.shf == base_result.shf
        @test result.lhf == base_result.lhf
        @test result.Cd == base_result.Cd

        # ustar derived from Cd matches
        @test isapprox(result.ustar, base_result.ustar; rtol = 1e-4)
    end

    @testset "Prescribed Fluxes and Cd with Defaults" begin
        # Test fallback logic: Pass nothing for T_sfc and q_sfc
        # Effectively T_sfc -> T_int, q_sfc -> q_int
        # The result won't match base_result exactly because T_sfc != T_int in base case.
        # But we verify it runs and produces reasonable output.

        specs = SF.FluxSpecs{FT}(
            shf = base_result.shf,
            lhf = base_result.lhf,
            Cd = base_result.Cd,
        )

        result = SF.surface_fluxes(
            param_set,
            T_int, q_int, FT(0), FT(0), ρ_int,
            nothing, nothing, # Pass nothing for guesses
            FT(0), Δz, FT(0),
            u_int, u_sfc,
            nothing,
            config,
            SF.PointValueScheme(),
            nothing,
            specs,
        )

        # Inputs pass through
        @test result.shf == base_result.shf
        @test result.lhf == base_result.lhf
        @test result.Cd == base_result.Cd

        # We can't strictly compare ustar/L_MO to base_result because ρ_sfc will differ 
        # (calculated from T_int instead of T_sfc).
        # Just check it returns valid numbers
        @test isfinite(result.ustar)
    end
end
