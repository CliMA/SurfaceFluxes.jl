using Test

import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD

const API_PRESSURE = 1.0e5

@inline function _density_from_state(thermo_params, T, pressure, qt)
    q_pt = TD.PhasePartition_equil_given_p(
        thermo_params,
        T,
        pressure,
        qt,
        TD.PhaseEquil,
    )
    return TD.air_density(thermo_params, T, pressure, q_pt)
end

@testset "Surface flux primitives API" begin
    FT = Float64
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
    thermo_params = SFP.thermodynamics_params(param_set)

    Tin = FT(300)
    qin = FT(0.01)
    Ts_guess = FT(298)
    qs_guess = FT(0.012)
    ρin = _density_from_state(thermo_params, Tin, FT(API_PRESSURE), qin)

    base_result = SF.surface_fluxes(
        param_set,
        Tin,
        qin,
        ρin,
        Ts_guess,
        qs_guess,
        FT(0),
        FT(10),
        FT(0),
    )

    @test base_result isa SF.SurfaceFluxConditions{FT}
    @test base_result.Cd > zero(FT)

    Ts_calls = Ref(0)
    qs_calls = Ref(0)
    update_Ts! = function (state, ctx)
        Ts_calls[] += 1
        return ctx.Tin - FT(5)
    end
    update_qs! = function (state, ctx)
        qs_calls[] += 1
        return ctx.qin + FT(2e-3)
    end

    hooked_result = SF.surface_fluxes(
        param_set,
        Tin,
        qin,
        ρin,
        Ts_guess,
        qs_guess,
        FT(0),
        FT(10),
        FT(0),
        nothing,
        nothing,
        nothing,
        SF.SurfaceFluxConfig(),
        SF.PointValueScheme(),
        nothing,
        nothing,
        update_Ts!,
        update_qs!,
    )

    @test Ts_calls[] > 0
    @test qs_calls[] > 0
    @test hooked_result.shf != base_result.shf

    config = SF.SurfaceFluxConfig(
        SF.charnock_momentum(; α = 0.02, scalar = 1e-4),
        SF.gustiness_constant(0.5),
    )

    config_result = SF.surface_fluxes(
        param_set,
        Tin,
        qin,
        ρin,
        Ts_guess,
        qs_guess,
        FT(0),
        FT(10),
        FT(0),
        nothing,
        nothing,
        nothing,
        config,
    )

    @test config_result.Cd != base_result.Cd
end
