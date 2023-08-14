using Test
using Logging

import QuadGK

import SurfaceFluxes
const SF = SurfaceFluxes
import SurfaceFluxes.UniversalFunctions as UF

import CLIMAParameters
const CP = CLIMAParameters

include(joinpath(pkgdir(SurfaceFluxes), "parameters", "create_parameters.jl"))
FT = Float64;
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
param_set = create_parameters(toml_dict, UF.BusingerType())
universal_functions(uft, L) = UF.universal_func(uft, L, create_uf_parameters(toml_dict, uft))

# TODO: Right now, we test these functions for
# type stability and correctness in the asymptotic
# limit. We may want to extend correctness tests.

@testset "UniversalFunctions" begin
    @testset "Type stability" begin
        FT = Float32
        ζ = FT(-2):FT(0.01):FT(200)
        for L in (-FT(10), FT(10))
            for uft in (
                UF.GryanikType(),
                UF.GrachevType(),
                UF.BusingerType(),
                UF.BeljaarsType(),
                UF.HoltslagType(),
                UF.ChengType(),
            )
                uf = universal_functions(uft, L)
                for transport in (UF.MomentumTransport(), UF.HeatTransport())
                    ϕ = UF.phi.(uf, ζ, transport)
                    @test eltype(ϕ) == FT
                    ψ = UF.psi.(uf, ζ, transport)
                    @test eltype(ψ) == FT
                end
            end
        end

        # More type stability (phi/psi):
        FT = Float32
        ζ = (-FT(1), FT(0.5) * eps(FT), 2 * eps(FT))
        for L in (-FT(10), FT(10))
            for uft in (
                UF.GryanikType(),
                UF.GrachevType(),
                UF.BusingerType(),
                UF.BeljaarsType(),
                UF.HoltslagType(),
                UF.ChengType(),
            )
                uf = universal_functions(uft, L)
                for transport in (UF.MomentumTransport(), UF.HeatTransport())
                    ϕ = UF.phi.(uf, ζ, transport)
                    @test eltype(ϕ) == FT
                    ψ = UF.psi.(uf, ζ, transport)
                    @test eltype(ψ) == FT
                end
            end
        end

        # More type stability (Psi):
        FT = Float32
        ζ = (-FT(1), -FT(0.5) * eps(FT), FT(0.5) * eps(FT), 2 * eps(FT))
        for L in (-FT(10), FT(10))
            for uft in (UF.GryanikType(), UF.BusingerType(), UF.BeljaarsType(), UF.HoltslagType(), UF.ChengType())
                uf = universal_functions(uft, L)
                for transport in (UF.MomentumTransport(), UF.HeatTransport())
                    Ψ = UF.Psi.(uf, ζ, transport)
                    @test eltype(Ψ) == FT
                end
            end
        end

        # Cheng formulation throws exception for Psi (ζ >> 0)
        ζ = FT(1):FT(0.01):FT(200)
        for L in (-FT(10), FT(10))
            uf = universal_functions(UF.ChengType(), L)
            for transport in (UF.MomentumTransport(), UF.HeatTransport())
                @test_logs (
                    :error,
                    "Volume-averaged form of Cheng (2005) SCF for ζ >> 0 not yet implemented. Use Nishizawa (2018) formulation instead.",
                ) UF.Psi.(uf, ζ, transport)
            end
        end

    end
    @testset "Asymptotic range" begin
        FT = Float32

        ϕ_h_ζ∞(uf::UF.Grachev, ζ) = 1 + FT(UF.b_h(uf))
        ϕ_m_ζ∞(uf::UF.Grachev, ζ) = FT(UF.a_m(uf)) / FT(UF.b_m(uf)) * ζ^FT(1 / 3)

        ϕ_h_ζ∞(uf::UF.Gryanik, ζ) = FT(1) + (FT(ζ) * FT(UF.Pr_0(uf)) * FT(UF.a_h(uf))) / (1 + FT(UF.b_h(uf)) * FT(ζ))
        ϕ_m_ζ∞(uf::UF.Gryanik, ζ) = FT(UF.a_m(uf) / UF.b_m(uf)^FT(2 / 3)) * ζ^FT(1 / 3)

        ϕ_h_ζ∞(uf::UF.Cheng, ζ) = 1 + FT(UF.a_h(uf))
        ϕ_m_ζ∞(uf::UF.Cheng, ζ) = 1 + FT(UF.a_m(uf))

        for L in (-FT(10), FT(10))
            for uft in (UF.GryanikType(), UF.GrachevType(), UF.ChengType())
                uf = universal_functions(uft, L)
                for ζ in FT(10) .^ (4, 6, 8, 10)
                    ϕ_h = UF.phi(uf, ζ, UF.HeatTransport())
                    @test isapprox(ϕ_h, ϕ_h_ζ∞(uf, ζ))
                end
                for ζ in FT(10) .^ (8, 9, 10)
                    ϕ_m = UF.phi(uf, ζ, UF.MomentumTransport())
                    @test isapprox(ϕ_m, ϕ_m_ζ∞(uf, ζ))
                end
            end
        end

    end

    # Test for Gryanik2021 Eq. 2 & 3; ensures Ψ(0) = 0
    @testset "Vanishes at Zero" begin
        FT = Float32
        for L in (-FT(10), FT(10))
            for uft in (UF.GryanikType(), UF.GrachevType())
                for transport in (UF.HeatTransport(), UF.MomentumTransport())
                    uf = universal_functions(uft, L)
                    Ψ_0 = UF.psi(uf, FT(0), transport)
                    @test isapprox(Ψ_0, FT(0))
                end
            end
        end
    end

    # Test that the integrated forms of the stability correction functions ψ(ζ) are consistent
    # when directly evaluated via numerical integration from ϕ(ζ).
    @testset "Test Correctness: ψ(ζ) = ∫ 𝒻(ϕ(ζ′)) dζ′" begin
        FloatType = (Float32, Float64)
        for FT in FloatType
            ζ_array = (FT(-20), FT(-10), FT(-1), -sqrt(eps(FT)), sqrt(eps(FT)), FT(1), FT(10), FT(20))
            for L in FT(10) .* sign.(ζ_array)
                for ζ in ζ_array
                    for uft in (
                        UF.GryanikType(),
                        UF.GrachevType(),
                        UF.BusingerType(),
                        UF.HoltslagType(),
                        UF.ChengType(),
                        UF.BeljaarsType(),
                    )
                        uf = universal_functions(uft, L)
                        for transport in (UF.MomentumTransport(), UF.HeatTransport())
                            # Compute ψ via numerical integration of 𝒻(ϕ(ζ))
                            ψ_int = QuadGK.quadgk(ζ′ -> (FT(1) - UF.phi(uf, ζ′, transport)) / ζ′, eps(FT), ζ)
                            # Compute ψ using function definitions of ψ(ζ)
                            ψ = UF.psi(uf, ζ, transport)
                            @test isapprox(ψ_int[1] - ψ, FT(0), atol = 10sqrt(eps(FT)))
                        end
                    end
                end
            end
        end
    end

end
