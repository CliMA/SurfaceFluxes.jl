# Tests type stability,continuity, and asymptotic behavior of universal functions.

using Test

import QuadGK
import SurfaceFluxes.UniversalFunctions as UF
import ClimaParams as CP

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

const TRANSPORTS = (UF.MomentumTransport(), UF.HeatTransport())

"""
    universal_parameter_sets(::Type{FT})

Return the available universal-function parameter structs for the requested
floating-point type. Each set is constructed using the default parameter values
from ClimaParams via the `UF.{Businger, Gryanik, Grachev}Params(FT)` constructors.
"""
universal_parameter_sets(::Type{FT}) where {FT <: AbstractFloat} =
    (UF.GryanikParams(FT), UF.GrachevParams(FT), UF.BusingerParams(FT))

# Some tests only apply to a subset of parameterizations (e.g., the integrated stability
# correction function `Psi` is not defined for the `Grachev` parameterization).
psi_parameter_sets(::Type{FT}) where {FT <: AbstractFloat} =
    (UF.GryanikParams(FT), UF.BusingerParams(FT))

# Asymptotic reference solutions for the similarity functions (nondimensional gradients)
ϕ_h_limit(p::UF.GrachevParams, ζ) = 1 + typeof(ζ)(p.b_h)
ϕ_m_limit(p::UF.GrachevParams, ζ) = typeof(ζ)(p.a_m / p.b_m) * ζ^(typeof(ζ)(1 / 3))

function ϕ_h_limit(p::UF.GryanikParams, ζ)
    Tζ = typeof(ζ)
    return Tζ(1) + (Tζ(ζ) * Tζ(p.Pr_0) * Tζ(p.a_h)) / (1 + Tζ(p.b_h) * Tζ(ζ))
end

ϕ_m_limit(p::UF.GryanikParams, ζ) =
    typeof(ζ)(p.a_m / p.b_m^(30 // 45)) * ζ^(typeof(ζ)(1 / 3)) # exponent = 2/3

"""
    neutral_velocity_increment(ufp, z, z0, L)

Dimensionless (κ / u★) velocity increment between heights `z0` and `z` using the
Monin–Obukhov similarity form. When `L → ∞` (neutral stratification) this must
collapse to the logarithmic law of the wall, i.e. `log(z / z0)`.
"""
function neutral_velocity_increment(ufp, z, z0, L)
    transport = UF.MomentumTransport()
    ζ = z / L
    ζ0 = z0 / L
    ψ_z = UF.psi(ufp, ζ, transport)
    ψ_ref = UF.psi(ufp, ζ0, transport)
    return log(z / z0) - ψ_z + ψ_ref
end

# ---------------------------------------------------------------------------
# Test-suite
# ---------------------------------------------------------------------------

@testset "UniversalFunctions" begin
    @testset "Type Stability (Phi & Psi)" begin
        for FT in (Float32, Float64)
            fine_grid = FT(-2):FT(0.01):FT(200)
            near_zero = (-FT(1), FT(0.5) * eps(FT), FT(2) * eps(FT))
            for ζ_values in (fine_grid, near_zero)
                for ufp in universal_parameter_sets(FT)
                    for transport in TRANSPORTS
                        ϕ = UF.phi.(Ref(ufp), ζ_values, Ref(transport))
                        ψ = UF.psi.(Ref(ufp), ζ_values, Ref(transport))
                        @test eltype(ϕ) == FT
                        @test eltype(ψ) == FT
                    end
                end
            end
        end
    end

    @testset "Type Stability (Psi)" begin
        for FT in (Float32, Float64)
            ζ_values = (-FT(1), -FT(0.5) * eps(FT), FT(0.5) * eps(FT), FT(2) * eps(FT))
            for ufp in psi_parameter_sets(FT)
                for transport in TRANSPORTS
                    Ψ = UF.Psi.(Ref(ufp), ζ_values, Ref(transport))
                    @test eltype(Ψ) == FT
                end
            end
        end
    end

    @testset "Neutral Logarithmic Velocity Profile" begin
        for FT in (Float32, Float64)
            z0 = FT(1)
            heights = (FT(2), FT(4), FT(8), FT(32), FT(128))
            # Use a very large Obukhov length to emulate neutral conditions.
            L_neutral = floatmax(FT) / FT(2)
            tol = max(FT(100) * eps(FT), FT(1e-12))
            for ufp in universal_parameter_sets(FT)
                for z in heights
                    Δu = neutral_velocity_increment(ufp, z, z0, L_neutral)
                    expected = log(z / z0)
                    @test isapprox(Δu, expected; atol = tol, rtol = FT(0))
                end
            end
        end
    end

    @testset "Asymptotic Behavior (|ζ| → ∞)" begin
        FT = Float32
        large_positive = FT(10) .^ (4, 6, 8, 10)
        very_large = FT(10) .^ (8, 9, 10)

        for ufp in (UF.GryanikParams(FT), UF.GrachevParams(FT))
            for ζ in large_positive
                ϕ_h = UF.phi(ufp, ζ, UF.HeatTransport())
                @test isapprox(ϕ_h, ϕ_h_limit(ufp, ζ); rtol = FT(5e-3))
            end
            for ζ in very_large
                ϕ_m = UF.phi(ufp, ζ, UF.MomentumTransport())
                @test isapprox(ϕ_m, ϕ_m_limit(ufp, ζ); rtol = FT(5e-3))
            end
        end
    end

    @testset "Monotonicity of ϕ(ζ)" begin
        # ϕ should increase with stability (ζ)
        for FT in (Float32, Float64)
            ζ_grid = range(FT(-5), FT(5), length = 100)

            for ufp in universal_parameter_sets(FT), transport in TRANSPORTS
                ϕ_vals = [UF.phi(ufp, ζ, transport) for ζ in ζ_grid]

                # Check that differences between consecutive points are positive
                diffs = diff(ϕ_vals)

                @test all(d -> d > -eps(FT), diffs) # Allow for numerical noise
            end
        end
    end

    @testset "ψ(0) == 0" begin
        FT = Float32
        for ufp in (UF.GryanikParams(FT), UF.GrachevParams(FT))
            for transport in TRANSPORTS
                ψ0 = UF.psi(ufp, FT(0), transport)
                @test isapprox(ψ0, FT(0); atol = eps(FT))
            end
        end
    end

    @testset "Neutral Continuity (ζ → 0)" begin
        for FT in (Float32, Float64)
            # Test close to zero to verify the limit, but not so close 
            # that we hit the linear approximation guards in the code.
            # 1e-6 is safe for Float32/64 to verify continuity behavior.
            ζ_small = FT(1e-6)

            # Tolerance must accommodate the slope times the step size.
            # Since we compare f(ε) to f(0), error is O(slope * ε).
            # Slope is ~5. 5 * 1e-6 = 5e-6.
            atol = FT(10 * ζ_small)

            for ufp in universal_parameter_sets(FT), transport in TRANSPORTS
                # 1. Limits for phi and psi (Available for ALL parameterizations)
                ϕ_0 = UF.phi(ufp, FT(0), transport)
                ψ_0 = UF.psi(ufp, FT(0), transport)

                # 2. Test Stable Branch (ζ > 0) -> Limit
                @test isapprox(UF.phi(ufp, ζ_small, transport), ϕ_0; atol = atol)
                @test isapprox(UF.psi(ufp, ζ_small, transport), ψ_0; atol = atol)

                # 3. Test Unstable Branch (ζ < 0) -> Limit
                @test isapprox(UF.phi(ufp, -ζ_small, transport), ϕ_0; atol = atol)
                @test isapprox(UF.psi(ufp, -ζ_small, transport), ψ_0; atol = atol)

                # 4. Limits for Psi (Available ONLY for Businger and Gryanik)
                # Grachev does not implement Psi, so we skip it to avoid MethodError
                if !(ufp isa UF.GrachevParams)
                    Ψ_0 = UF.Psi(ufp, FT(0), transport)
                    @test isapprox(UF.Psi(ufp, ζ_small, transport), Ψ_0; atol = atol)
                    @test isapprox(UF.Psi(ufp, -ζ_small, transport), Ψ_0; atol = atol)
                end
            end

        end
    end

    @testset "Derivative Consistency ϕ(ζ) ≈ ϕ(0) - ζ·ψ'(ζ)" begin
        # Test that the analytical psi is consistent with phi via finite differences
        for FT in (Float32, Float64)
            ζ_samples = (FT(-5), FT(-1), FT(-0.1), FT(0.1), FT(1), FT(5))
            ϵ = cbrt(eps(FT)) # Finite difference step size

            for ζ in ζ_samples, ufp in universal_parameter_sets(FT), transport in TRANSPORTS
                # Central difference approximation of ψ'(ζ)
                ψ_plus = UF.psi(ufp, ζ + ϵ, transport)
                ψ_minus = UF.psi(ufp, ζ - ϵ, transport)
                dψ_dζ = (ψ_plus - ψ_minus) / (2ϵ)

                ϕ_0 = UF.phi(ufp, FT(0), transport)
                ϕ_expected = ϕ_0 - ζ * dψ_dζ
                ϕ_actual = UF.phi(ufp, ζ, transport)

                # Tolerances need to be loose for FD, especially with changing curvature
                @test isapprox(ϕ_actual, ϕ_expected; rtol = sqrt(ϵ))
            end
        end
    end

    @testset "Integral Consistency ψ(ζ) = ∫(ϕ(0)-ϕ)/ζ′ dζ′" begin
        for FT in (Float32, Float64)
            ζ_samples = (
                FT(-20),
                FT(-10),
                FT(-1),
                -sqrt(eps(FT)),
                sqrt(eps(FT)),
                FT(1),
                FT(10),
                FT(20),
            )
            for ζ in ζ_samples, ufp in universal_parameter_sets(FT), transport in TRANSPORTS
                # Determine the neutral limit ϕ(0) dynamically
                # For Businger this is 1. For Gryanik/Grachev heat this is Pr_0.
                ϕ_0 = UF.phi(ufp, FT(0), transport)

                # quadgk returns (value, error)
                ψ_int, _ = QuadGK.quadgk(
                    ζ′ -> (ϕ_0 - UF.phi(ufp, ζ′, transport)) / ζ′,
                    eps(FT),
                    ζ;
                    rtol = 1e-9,
                )
                ψ = UF.psi(ufp, ζ, transport)
                @test isapprox(ψ_int, ψ; atol = FT(10) * sqrt(eps(FT)))
            end
        end
    end

    @testset "Small-ζ Linearization Consistency" begin
        # Checks that the implementation of Psi near zero (using linear approximations)
        # is consistent with the theoretical slope derived from phi parameters.
        # The limit of Psi(ζ)/ζ as ζ->0 should be -phi'(0)/2.
        # For most functions here: phi(ζ) ~ phi(0) + a*ζ  =>  Psi(ζ) ~ -a*ζ/2

        for FT in (Float32, Float64)
            # Choose a ζ small enough to trigger the `abs(ζ) < eps(FT)` guard
            # inside the Psi functions.
            ζ_tiny = FT(0.9) * eps(FT)

            for ufp in psi_parameter_sets(FT) # Skip Grachev as Psi not implemented

                # 1. Momentum
                # phi_m ~ 1 + a_m * ζ
                # Psi_m ~ -a_m * ζ / 2
                Psi_m = UF.Psi(ufp, ζ_tiny, UF.MomentumTransport())
                expected_m = -FT(ufp.a_m) * ζ_tiny / 2
                @test isapprox(Psi_m, expected_m; rtol = sqrt(eps(FT)))

                # 2. Heat
                # phi_h ~ Pr_0 * (1 + a_h * ζ)  [Gryanik] or 1 + a_h * ζ [Businger]
                # The linear coefficient 'a' depends on the parameterization structure.
                Psi_h = UF.Psi(ufp, ζ_tiny, UF.HeatTransport())

                if ufp isa UF.GryanikParams
                    # Gryanik Heat: phi_h ~ Pr_0 + Pr_0 * a_h * ζ
                    # Slope is Pr_0 * a_h
                    expected_h = -FT(ufp.Pr_0) * FT(ufp.a_h) * ζ_tiny / 2
                else
                    # Businger Heat: phi_h ~ Pr_0 + a_h * ζ
                    # Psi_h ~ -a_h * ζ / 2
                    expected_h = -FT(ufp.a_h) * ζ_tiny / 2
                end
                @test isapprox(Psi_h, expected_h; rtol = sqrt(eps(FT)))
            end
        end
    end

    @testset "Bulk Richardson Number" begin
        for FT in (Float32, Float64)
            # Choose a grid that avoids exactly 0 for monotonicity check steps if needed, 
            # though we test 0 explicitly.
            # Range including stable and unstable
            ζ_grid = range(FT(-10), FT(10), length = 200)
            Δz = FT(10)
            z0m = FT(0.1)
            z0h = FT(0.01)

            schemes = (UF.PointValueScheme(), UF.LayerAverageScheme())

            for ufp in universal_parameter_sets(FT), scheme in schemes
                # Skip Grachev for LayerAverageScheme as Psi is not implemented
                if scheme isa UF.LayerAverageScheme && ufp isa UF.GrachevParams
                    continue
                end

                # 1. Neutral limit: Ri_b(0) should be 0
                # Because thermal stratification is zero, buoyancy production is zero.
                @test isapprox(
                    UF.bulk_richardson_number(ufp, Δz, FT(0), z0m, z0h, scheme),
                    FT(0);
                    atol = eps(FT),
                )

                # 2. Continuity near neutral limit
                ε = sqrt(eps(FT))
                Rib_pos = UF.bulk_richardson_number(ufp, Δz, ε, z0m, z0h, scheme)
                Rib_neg = UF.bulk_richardson_number(ufp, Δz, -ε, z0m, z0h, scheme)
                # Should be small and order of ε
                @test isapprox(Rib_pos, FT(0); atol = 10ε)
                @test isapprox(Rib_neg, FT(0); atol = 10ε)

                # 3. Monotonicity in ζ
                # Ri_b should generally increase with ζ.
                # Ri_b ~ ζ * F_h / F_m^2. 
                # F_h and F_m are positive and monotonic.
                # We check if Ri_b(ζ_{i+1}) > Ri_b(ζ_i).

                # Compute Ri_b across the grid
                Ris = [
                    UF.bulk_richardson_number(ufp, Δz, ζ, z0m, z0h, scheme) for ζ in ζ_grid
                ]

                # Check sorted
                @test issorted(Ris)
            end
        end
    end

    @testset "Dimensionless Profile" begin
        for FT in (Float32, Float64)
            z0 = FT(0.1)
            Δz = FT(10)
            ζ_grid = range(FT(-2), FT(2), length = 20)

            schemes = (UF.PointValueScheme(), UF.LayerAverageScheme())

            for ufp in universal_parameter_sets(FT), transport in TRANSPORTS,
                scheme in schemes
                # Skip Grachev for LayerAverageScheme
                if scheme isa UF.LayerAverageScheme && ufp isa UF.GrachevParams
                    continue
                end

                # 1. Neutral limit (ζ -> 0)
                F_neutral = UF.dimensionless_profile(ufp, Δz, FT(0), z0, transport, scheme)
                slope = UF.phi(ufp, FT(0), transport)

                expected_neutral = if scheme isa UF.PointValueScheme
                    slope * log(Δz / z0)
                else # LayerAverageScheme
                    # N&K 2018 approximation: slope * (log(Δz/z0) - (1 - z0/Δz))
                    slope * (log(Δz / z0) - FT(1) + z0 / Δz)
                end
                @test isapprox(F_neutral, expected_neutral; atol = 10 * eps(FT))

                # 2. Continuity near neutral limit
                ε = sqrt(eps(FT))
                F_pos = UF.dimensionless_profile(ufp, Δz, ε, z0, transport, scheme)
                F_neg = UF.dimensionless_profile(ufp, Δz, -ε, z0, transport, scheme)
                @test isapprox(F_pos, F_neutral; atol = 10ε)
                @test isapprox(F_neg, F_neutral; atol = 10ε)

                # 3. Monotonicity in ζ
                # Dimensionless profile F(ζ) should generally increase with ζ (more stable = larger gradient)
                Fs = [
                    UF.dimensionless_profile(ufp, Δz, ζ, z0, transport, scheme) for
                    ζ in ζ_grid
                ]
                @test issorted(Fs)
            end
        end
    end
end
