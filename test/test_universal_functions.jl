using Test

import QuadGK

import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import ClimaParams as CP
import SurfaceFluxes.Parameters as SFP

FT = Float32
param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
thermo_params = param_set.thermo_params

# TODO: Right now, we test these functions for
# type stability and correctness in the asymptotic
# limit. We may want to extend correctness tests.

@testset "UniversalFunctions" begin
    @testset "Type stability" begin
        FT = Float32
        ζ = FT(-2):FT(0.01):FT(200)
        for ufp in (
            UF.GryanikParams(FT),
            UF.GrachevParams(FT),
            UF.BusingerParams(FT),
        )
            for transport in (UF.MomentumTransport(), UF.HeatTransport())
                ϕ = UF.phi.(ufp, ζ, transport)
                @test eltype(ϕ) == FT
                ψ = UF.psi.(ufp, ζ, transport)
                @test eltype(ψ) == FT
            end
        end

        # More type stability (phi/psi):
        FT = Float32
        ζ = (-FT(1), FT(0.5) * eps(FT), 2 * eps(FT))
        for ufp in (
            UF.GryanikParams(FT),
            UF.GrachevParams(FT),
            UF.BusingerParams(FT),
        )
            for transport in (UF.MomentumTransport(), UF.HeatTransport())
                ϕ = UF.phi.(ufp, ζ, transport)
                @test eltype(ϕ) == FT
                ψ = UF.psi.(ufp, ζ, transport)
                @test eltype(ψ) == FT
            end
        end

        # More type stability (Psi):
        FT = Float32
        ζ = (-FT(1), -FT(0.5) * eps(FT), FT(0.5) * eps(FT), 2 * eps(FT))
        for ufp in (UF.GryanikParams(FT), UF.BusingerParams(FT))
            for transport in (UF.MomentumTransport(), UF.HeatTransport())
                Ψ = UF.Psi.(ufp, ζ, transport)
                @test eltype(Ψ) == FT
            end
        end

    end
    @testset "Asymptotic range" begin
        FT = Float32

        ϕ_h_ζ∞(p::UF.GrachevParams, ζ) = 1 + FT(p.b_h)
        ϕ_m_ζ∞(p::UF.GrachevParams, ζ) = FT(p.a_m) / FT(p.b_m) * ζ^FT(1 / 3)

        ϕ_h_ζ∞(p::UF.GryanikParams, ζ) =
            FT(1) + (FT(ζ) * FT(p.Pr_0) * FT(p.a_h)) / (1 + FT(p.b_h) * FT(ζ))
        ϕ_m_ζ∞(p::UF.GryanikParams, ζ) =
            FT(p.a_m / p.b_m^FT(2 / 3)) * ζ^FT(1 / 3)


        for ufp in (UF.GryanikParams(FT), UF.GrachevParams(FT))
            for ζ in FT(10) .^ (4, 6, 8, 10)
                ϕ_h = UF.phi(ufp, ζ, UF.HeatTransport())
                @test isapprox(ϕ_h, ϕ_h_ζ∞(ufp, ζ))
            end
            for ζ in FT(10) .^ (8, 9, 10)
                ϕ_m = UF.phi(ufp, ζ, UF.MomentumTransport())
                @test isapprox(ϕ_m, ϕ_m_ζ∞(ufp, ζ))
            end
        end

    end

    # Test for Gryanik2021 Eq. 2 & 3; ensures Ψ(0) = 0
    @testset "Vanishes at Zero" begin
        FT = Float32
        for ufp in (UF.GryanikParams(FT), UF.GrachevParams(FT))
            for transport in (UF.HeatTransport(), UF.MomentumTransport())
                Ψ_0 = UF.psi(ufp, FT(0), transport)
                @test isapprox(Ψ_0, FT(0), atol = eps(FT))
            end
        end
    end

    # Test that the integrated forms of the stability correction functions ψ(ζ) are consistent
    # when directly evaluated via numerical integration from ϕ(ζ).
    @testset "Test Correctness: ψ(ζ) = ∫ 𝒻(ϕ(ζ′)) dζ′" begin
        FloatType = (Float32, Float64)
        for FT in FloatType
            ζ_array = (
                FT(-20),
                FT(-10),
                FT(-1),
                -sqrt(eps(FT)),
                sqrt(eps(FT)),
                FT(1),
                FT(10),
                FT(20),
            )
            for ζ in ζ_array
                for ufp in (
                    UF.GryanikParams(FT),
                    UF.GrachevParams(FT),
                    UF.BusingerParams(FT),
                )
                    for transport in (UF.MomentumTransport(), UF.HeatTransport())
                        ψ_int = QuadGK.quadgk(
                            ζ′ -> (FT(1) - UF.phi(ufp, ζ′, transport)) / ζ′,
                            eps(FT),
                            ζ,
                        )
                        ψ = UF.psi(ufp, ζ, transport)
                        @test isapprox(
                            ψ_int[1] - ψ,
                            FT(0),
                            atol = 10sqrt(eps(FT)),
                        )
                    end
                end
            end
        end
    end

end
