using Test

import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import ClimaParams as CP
import Thermodynamics as TD
using RootSolvers
const RS = RootSolvers

# Reuse helper functions from test_convergence.jl
include("test_convergence.jl")

@testset "RootSolvers.jl integration tests" begin
    # Test cases: stable and unstable conditions
    test_cases = [
        # (T_sfc, T_in, RH_in, Δqt, z, speed, wind_dir, z0m, z0b_factor, p, description)
        (282.0, 280.0, 0.85, 0.0, 20.0, 5.0, (1.0, 0.0), 1e-3, 1.0, 1.0e5, "stable"),
        (280.0, 282.0, 0.85, 0.0, 20.0, 5.0, (1.0, 0.0), 1e-3, 1.0, 1.0e5, "unstable"),
        (282.0, 282.0, 0.85, 0.0, 20.0, 5.0, (1.0, 0.0), 1e-3, 1.0, 1.0e5, "neutral"),
    ]

    for FT in (Float32, Float64)
        param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
        scheme = SF.PointValueScheme()
        tol = FT(1e-6)
        maxiter = 20

        for (T_sfc, T_in, RH_in, Δqt, z, speed, wind_dir, z0m, z0b_factor, p, desc) in test_cases
            # Build surface condition
            qt_in = qt_from_RH(FT, RH_in, T_in, p)
            qt_sfc = FT(clamp(qt_in + Δqt, 0.0, 1.0))
            wind = (FT(speed * wind_dir[1]), FT(speed * wind_dir[2]))
            z0b = FT(max(z0m * z0b_factor, 1e-6))

            case = (
                T_sfc = FT(T_sfc),
                T_in = FT(T_in),
                qt_in = qt_in,
                qt_sfc = qt_sfc,
                z = FT(z),
                wind = wind,
                z0m = FT(z0m),
                z0b = z0b,
                pressure = FT(p),
            )

            sc = build_surface_condition(param_set, case, SF.ScalarRoughness())

            # Create initial guess for similarity scales
            δ = sign(SF.ΔDSEᵥ(param_set, sc))
            u★₀ = FT(0.1)
            𝓁u₀ = SF.compute_z0(
                u★₀,
                param_set,
                sc,
                sc.roughness_model,
                UF.MomentumTransport(),
            )
            𝓁θ₀ = SF.compute_z0(
                u★₀,
                param_set,
                sc,
                sc.roughness_model,
                UF.HeatTransport(),
            )
            𝓁q₀ = SF.compute_z0(
                u★₀,
                param_set,
                sc,
                sc.roughness_model,
                UF.HeatTransport(),
            )

            if SF.ΔDSEᵥ(param_set, sc) >= FT(0)
                X★₀ = (
                    u★ = u★₀,
                    DSEᵥ★ = FT(δ),
                    θᵥ★ = FT(δ),
                    q★ = FT(δ),
                    L★ = FT(10),
                    𝓁u = 𝓁u₀,
                    𝓁θ = 𝓁θ₀,
                    𝓁q = 𝓁q₀,
                )
            else
                X★₀ = (
                    u★ = u★₀,
                    DSEᵥ★ = FT(δ),
                    θᵥ★ = FT(δ),
                    q★ = FT(δ),
                    L★ = FT(-10),
                    𝓁u = 𝓁u₀,
                    𝓁θ = 𝓁θ₀,
                    𝓁q = 𝓁q₀,
                )
            end

            # Helper function to track iterations for fixed point method
            function fixed_point_with_tracking()
                X★_prev = X★₀
                qₛ = SF.surface_specific_humidity(param_set, sc)
                iter_count = 0
                
                function fixed_point_func(X★_in)
                    return SF.iterate_interface_fluxes(
                        sc,
                        qₛ,
                        X★_in,
                        SF.ts_in(sc),
                        SF.ts_sfc(sc),
                        scheme,
                        param_set,
                    )
                end
                
                function convergence_check(X★_prev, X★_curr)
                    Ri_b_prev = SF.bulk_richardson_number(param_set, sc, X★_prev, scheme)
                    Ri_b_curr = SF.bulk_richardson_number(param_set, sc, X★_curr, scheme)
                    return abs(Ri_b_curr - Ri_b_prev) ≤ tol
                end
                
                for iter in 1:maxiter
                    iter_count = iter
                    X★_curr = fixed_point_func(X★_prev)
                    if convergence_check(X★_prev, X★_curr)
                        return (X★_curr, iter_count)
                    end
                    X★_prev = X★_curr
                end
                return (X★_prev, iter_count)
            end

            # Test FixedPointIteration (baseline) with timing
            elapsed_fixed = @elapsed begin
                result_fixed = fixed_point_with_tracking()
            end
            X★_fixed, iter_fixed = result_fixed

            # Compute brackets for RootSolvers methods based on initial guess
            # This matches the logic in obukhov_iteration_rootsolver
            ζ₀ = SF.Δz(sc) / X★₀.L★
            if ζ₀ >= 0
                # Stable conditions: ζ > 0
                bracket_low = max(FT(1e-6), ζ₀ * FT(0.1))
                bracket_high = min(FT(10.0), max(FT(1.0), ζ₀ * FT(10.0)))
            else
                # Unstable conditions: ζ < 0
                bracket_low = max(FT(-10.0), min(FT(-1.0), ζ₀ * FT(10.0)))
                bracket_high = min(FT(-1e-6), ζ₀ * FT(0.1))
            end

            # Helper to get iteration count from RootSolvers solution
            # We need to access the solution object, so we'll create a wrapper
            function rootsolver_with_tracking(method)
                # RootSolvers doesn't expose iteration count easily in the public API
                # We'll use maxiter as an upper bound and note that it converged
                sol = (converged = true, iters = maxiter)  # Placeholder
                X★_result = nothing
                elapsed = @elapsed begin
                    # Call obukhov_iteration and measure time
                    X★_result = SF.obukhov_iteration(
                        X★₀,
                        sc,
                        scheme,
                        param_set,
                        tol,
                        maxiter,
                        method,
                    )
                end
                return (X★_result, sol, elapsed)
            end

            # Test BrentsMethod with RootSolvers.BrentsMethod
            X★_brent, sol_brent, elapsed_brent = rootsolver_with_tracking(
                RS.BrentsMethod(bracket_low, bracket_high),
            )
            iter_brent = sol_brent.iters  # Will be maxiter as placeholder

            # Test SecantMethod with RootSolvers.SecantMethod
            X★_secant, sol_secant, elapsed_secant = rootsolver_with_tracking(
                RS.SecantMethod(bracket_low, bracket_high),
            )
            iter_secant = sol_secant.iters  # Will be maxiter as placeholder

            # Verify all methods produce valid results
            @testset "$desc conditions (FT=$FT)" begin
                # Check that all methods produce finite, reasonable values
                for (method_name, X★) in [
                    ("FixedPoint", X★_fixed),
                    ("Brent", X★_brent),
                    ("Secant", X★_secant),
                ]
                    @testset "$method_name: validity checks" begin
                        @test isfinite(X★.u★)
                        @test X★.u★ > FT(0)
                        @test isfinite(X★.L★)
                        @test X★.L★ != FT(0)
                        @test isfinite(X★.DSEᵥ★)
                        @test isfinite(X★.θᵥ★)
                        @test isfinite(X★.q★)
                        @test X★.𝓁u > FT(0)
                        @test X★.𝓁θ > FT(0)
                        @test X★.𝓁q > FT(0)
                    end
                end

                # Check that all methods converge to similar Richardson numbers
                Ri_b_fixed = SF.bulk_richardson_number(param_set, sc, X★_fixed, scheme)
                Ri_b_brent = SF.bulk_richardson_number(param_set, sc, X★_brent, scheme)
                Ri_b_secant = SF.bulk_richardson_number(param_set, sc, X★_secant, scheme)

                # All methods should produce similar Richardson numbers
                # (within tolerance, accounting for different convergence paths)
                rtol = FT(1e-3)  # Relative tolerance for comparison
                @testset "Richardson number consistency" begin
                    @test isapprox(
                        Ri_b_brent,
                        Ri_b_fixed,
                        rtol = rtol,
                    )
                    @test isapprox(
                        Ri_b_secant,
                        Ri_b_fixed,
                        rtol = rtol,
                    )
                end

                # Check that stability parameter ζ is consistent
                ζ_fixed = SF.Δz(sc) / X★_fixed.L★
                ζ_brent = SF.Δz(sc) / X★_brent.L★
                ζ_secant = SF.Δz(sc) / X★_secant.L★

                @testset "Stability parameter consistency" begin
                    @test isapprox(ζ_brent, ζ_fixed, rtol = rtol)
                    @test isapprox(ζ_secant, ζ_fixed, rtol = rtol)
                end

                # Check that friction velocity is similar
                @testset "Friction velocity consistency" begin
                    @test isapprox(
                        X★_brent.u★,
                        X★_fixed.u★,
                        rtol = rtol,
                    )
                    @test isapprox(
                        X★_secant.u★,
                        X★_fixed.u★,
                        rtol = rtol,
                    )
                end

                # Performance summary
                @info "Performance summary for $desc conditions (FT=$FT):" iterations = (
                    FixedPoint = iter_fixed,
                    Brent = iter_brent,
                    Secant = iter_secant,
                ) time_seconds = (
                    FixedPoint = elapsed_fixed,
                    Brent = elapsed_brent,
                    Secant = elapsed_secant,
                )
            end
        end
    end
end

@testset "RootSolvers method type stability" begin
    # Test that FixedPointIteration is properly typed
    @test SF.FixedPointIteration() isa SF.SolverMethod

    # Test that RootSolvers methods are properly typed
    FT = Float64
    bracket_low = FT(1e-6)
    bracket_high = FT(10.0)
    @test RS.BrentsMethod(bracket_low, bracket_high) isa RS.BrentsMethod
    @test RS.SecantMethod(bracket_low, bracket_high) isa RS.SecantMethod

    # Test that RootSolvers methods can be passed to obukhov_iteration
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
    scheme = SF.PointValueScheme()

    # Create a simple test case
    case = (
        T_sfc = FT(282.0),
        T_in = FT(280.0),
        qt_in = FT(0.01),
        qt_sfc = FT(0.01),
        z = FT(20.0),
        wind = (FT(5.0), FT(0.0)),
        z0m = FT(1e-3),
        z0b = FT(1e-3),
        pressure = FT(1.0e5),
    )
    sc = build_surface_condition(param_set, case, SF.ScalarRoughness())

    δ = sign(SF.ΔDSEᵥ(param_set, sc))
    u★₀ = FT(0.1)
    𝓁u₀ = SF.compute_z0(
        u★₀,
        param_set,
        sc,
        sc.roughness_model,
        UF.MomentumTransport(),
    )
    𝓁θ₀ = SF.compute_z0(
        u★₀,
        param_set,
        sc,
        sc.roughness_model,
        UF.HeatTransport(),
    )
    𝓁q₀ = SF.compute_z0(
        u★₀,
        param_set,
        sc,
        sc.roughness_model,
        UF.HeatTransport(),
    )

    X★₀ = (
        u★ = u★₀,
        DSEᵥ★ = FT(δ),
        θᵥ★ = FT(δ),
        q★ = FT(δ),
        L★ = FT(10),
        𝓁u = 𝓁u₀,
        𝓁θ = 𝓁θ₀,
        𝓁q = 𝓁q₀,
    )

    tol = FT(1e-6)
    maxiter = 20

    # Compute brackets for RootSolvers methods
    ζ₀ = SF.Δz(sc) / X★₀.L★
    if ζ₀ >= 0
        bracket_low = max(FT(1e-6), ζ₀ * FT(0.1))
        bracket_high = min(FT(10.0), max(FT(1.0), ζ₀ * FT(10.0)))
    else
        bracket_low = max(FT(-10.0), min(FT(-1.0), ζ₀ * FT(10.0)))
        bracket_high = min(FT(-1e-6), ζ₀ * FT(0.1))
    end

    # Test that each method type works and returns NamedTuple
    result_fixed = SF.obukhov_iteration(
        X★₀,
        sc,
        scheme,
        param_set,
        tol,
        maxiter,
        SF.FixedPointIteration(),
    )
    @test result_fixed isa NamedTuple
    @test haskey(result_fixed, :u★)
    @test haskey(result_fixed, :L★)

    # Test RootSolvers.BrentsMethod
    result_brent = SF.obukhov_iteration(
        X★₀,
        sc,
        scheme,
        param_set,
        tol,
        maxiter,
        RS.BrentsMethod(bracket_low, bracket_high),
    )
    @test result_brent isa NamedTuple
    @test haskey(result_brent, :u★)
    @test haskey(result_brent, :L★)

    # Test RootSolvers.SecantMethod
    result_secant = SF.obukhov_iteration(
        X★₀,
        sc,
        scheme,
        param_set,
        tol,
        maxiter,
        RS.SecantMethod(bracket_low, bracket_high),
    )
    @test result_secant isa NamedTuple
    @test haskey(result_secant, :u★)
    @test haskey(result_secant, :L★)

end

