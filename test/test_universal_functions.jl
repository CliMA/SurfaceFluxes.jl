using Test

import SurfaceFluxes
const SF = SurfaceFluxes
const UF = SF.UniversalFunctions


# create the parameter sets for Gryanik and Grachev from src_parameter_dict
import CLIMAParameters
import Thermodynamics.ThermodynamicsParameters
import SurfaceFluxes.SurfaceFluxesParameters
import SurfaceFluxes.GryanikParameters
import SurfaceFluxes.GrachevParameters

gryanik_param_set = SurfaceFluxesParameters(
    src_parameter_dict,
    Gryanik(src_parameter_dict)
    ThermodynamicsParameters(src_parameter_dict),
)

grachev_param_set = SurfaceFluxesParameters(
    src_parameter_dict,
    Grachev(src_parameter_dict)
    ThermodynamicsParameters(src_parameter_dict),
)


# TODO: Right now, we test these functions for
# type stability and correctness in the asymptotic
# limit. We may want to extend correctness tests.

@testset "UniversalFunctions" begin
    @testset "Type stability" begin
        FT = Float32
        ζ = FT(-2):FT(0.01):FT(200)
        for L in (-FT(10), FT(10))
            for uf in (UF.Gryanik(L, gryanik_param_set), UF.Grachev(L, grachev_param_set), UF.Businger(L, businger_param_set))
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
            for uf in (UF.Gryanik(L, gryanik_param_set), UF.Grachev(L, grachev_param_set), UF.Businger(L, businger_param_set))
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
            uf = UF.Businger(L,businger_param_set)
            for transport in (UF.MomentumTransport(), UF.HeatTransport())
                Ψ = UF.Psi.(uf, ζ, transport)
                @test eltype(Ψ) == FT
            end
        end
    end
    @testset "Conversions" begin
        FT = Float32
        ζ = FT(10)
        L = FT(10)
        
        uf = UF.Gryanik(L, gryanik_param_set)
        @test UF.Businger(uf) isa UF.Businger
        @test UF.Grachev(uf) isa UF.Grachev

        uf = UF.Grachev(L, grachev_param_set)
        @test UF.Businger(uf) isa UF.Businger
        @test UF.Gryanik(uf) isa UF.Gryanik

        uf = UF.Businger(L, businger_param_set)
        @test UF.Grachev(uf) isa UF.Grachev
        @test UF.Gryanik(uf) isa UF.Gryanik
    end
    @testset "Asymptotic range" begin
        FT = Float32

        ϕ_h_ζ∞(uf::UF.Grachev) = 1 + FT(uf.param_set.b_h)
        ϕ_m_ζ∞(uf::UF.Grachev, ζ) = FT(uf.param_set.a_m) / FT(uf.param_set.b_m) * ζ^FT(1 / 3)

        ϕ_h_ζ∞(uf::UF.Gryanik) = FT(uf.param_set.Pr_0) * (1 + FT(uf.param_set.a_h / uf.param_set.b_h))
        ϕ_m_ζ∞(uf::UF.Gryanik, ζ) = FT(uf.param_set.a_m / uf.param_set.b_m^FT(2 / 3)) * ζ^FT(1 / 3)

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
