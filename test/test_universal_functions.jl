using Test

import SurfaceFluxes
const SF = SurfaceFluxes
const UF = SF.UniversalFunctions

import CLIMAParameters
const CP = CLIMAParameters

struct EarthParameterSet <: CP.AbstractEarthParameterSet end
const param_set = EarthParameterSet()

# TODO: Right now, we test these functions for
# type stability and correctness in the asymptotic
# limit. We may want to extend correctness tests.

@testset "UniversalFunctions" begin
    @testset "Type stability" begin
        FT = Float32
        ζ = FT(-2):FT(0.01):FT(200)
        for L in (-FT(10), FT(10))
            args = (L, param_set)
            for uf in (UF.Gryanik(args...), UF.Grachev(args...), UF.Businger(args...))
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
            args = (L, param_set)
            for uf in (UF.Gryanik(args...), UF.Grachev(args...), UF.Businger(args...))
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
            args = (L, param_set)
            uf = UF.Businger(args...)
            for transport in (UF.MomentumTransport(), UF.HeatTransport())
                Ψ = UF.Psi.(uf, ζ, transport)
                @test eltype(Ψ) == FT
            end
        end
    end
    @testset "Asymptotic range" begin
        FT = Float32

        ϕ_h_ζ∞(uf::UF.Grachev) = 1 + FT(UF.b_h(uf))
        ϕ_m_ζ∞(uf::UF.Grachev, ζ) = FT(UF.a_m(uf)) / FT(UF.b_m(uf)) * ζ^FT(1 / 3)

        ϕ_h_ζ∞(uf::UF.Gryanik) = FT(UF.Pr_0(uf)) * (1 + FT(UF.a_h(uf) / UF.b_h(uf)))
        ϕ_m_ζ∞(uf::UF.Gryanik, ζ) = FT(UF.a_m(uf) / UF.b_m(uf)^FT(2 / 3)) * ζ^FT(1 / 3)

        for L in (-FT(10), FT(10))
            args = (L, param_set)
            for uf in (UF.Grachev(args...), UF.Gryanik(args...))
                for ζ in FT(10) .^ (4, 6, 8, 10)
                    ϕ_h = UF.phi(uf, ζ, UF.HeatTransport())
                    @test isapprox(ϕ_h, ϕ_h_ζ∞(uf))
                end
                for ζ in FT(10) .^ (8, 9, 10)
                    ϕ_m = UF.phi(uf, ζ, UF.MomentumTransport())
                    @test isapprox(ϕ_m, ϕ_m_ζ∞(uf, ζ))
                end
            end
        end

    end

end
