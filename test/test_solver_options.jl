using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import ClimaParams as CP

# We need to construct a case where convergence usually takes a few iterations,
# then verify that using forced_fixed_iters=true works (doesn't crash) and produces similar results.
# Note: Verifying exact iteration count requires internal instrumentation,
# but we can verify that it runs and converges to a reasonable value.

@testset "SurfaceFluxes Forced Fixed Iters" begin
    FT = Float64
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
    thermo_params = SFP.thermodynamics_params(param_set)

    # Define a standard case
    T_sfc = FT(300)
    ρ_sfc = FT(1.2)
    T_int = FT(295) # Unstable
    ρ_int = FT(1.1)
    u_int = (FT(5), FT(0))
    u_sfc = (FT(0), FT(0))
    z_int = FT(10)
    z0m = FT(0.01)
    z0h = FT(0.01)

    inputs = SF.build_surface_flux_inputs(
        T_int, FT(0), FT(0), FT(0), ρ_int,
        T_sfc, FT(0), FT(0), z_int, FT(0),
        u_int, u_sfc,
        SF.SurfaceFluxConfig(
            SF.ConstantRoughnessParams(FT(z0m), FT(z0h)),
            SF.ConstantGustinessSpec(FT(1.0)),
            SF.DryModel(),
        ),
        nothing,
        SF.FluxSpecs{FT}(),
        nothing, nothing,
    )

    # 1. Standard run
    opts_std = SF.SolverOptions{FT}(maxiter = 10, tol = 1e-3, forced_fixed_iters = false)
    res_std = SF.solve_monin_obukhov(param_set, inputs, SF.PointValueScheme(), opts_std)

    # 2. Fixed iter run
    opts_fixed = SF.SolverOptions{FT}(maxiter = 10, tol = 1e-3, forced_fixed_iters = true)
    res_fixed = SF.solve_monin_obukhov(param_set, inputs, SF.PointValueScheme(), opts_fixed)

    # They should be very close (since standard usually converges in <10 iters and fixed continues to refine or holds)
    @test isapprox(res_std.shf, res_fixed.shf; atol = 1e-2)
    @test isapprox(res_std.ustar, res_fixed.ustar; atol = 1e-3)
end
