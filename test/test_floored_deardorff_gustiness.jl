# Floored, closed-form Deardorff gustiness: at fixed ζ the closed form is a fixed point
# of the bulk relations for both discretization schemes, and the floor applies when the
# surface is not warmer than the air. Through the full solve, the returned friction
# velocity is that of the post-solve Deardorff gustiness, the free-convection limit (zero
# wind) is well posed, stable conditions reduce to the constant floor, and Float32 and
# dual numbers propagate.

module TestFlooredDeardorffGustiness

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import ClimaParams as CP
import ForwardDiff

function make_inputs(FT, param_set, gust; U = 2, ΔT = 10, Δz = 10, z0m = 0.01, q = 0.01)
    return SF.build_surface_flux_inputs(
        FT(300), FT(q), FT(0), FT(0), FT(1.2), FT(300 + ΔT), FT(q), FT(0), FT(Δz), FT(0),
        (FT(U), FT(0)), (FT(0), FT(0)),
        SF.SurfaceFluxConfig(SF.ConstantRoughnessParams(FT(z0m), FT(z0m / 10)), gust),
        nothing, SF.FluxSpecs{FT}(), nothing, nothing,
    )
end

@testset "FlooredDeardorffGustinessSpec" begin
    FT = Float64
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
    β = SFP.gustiness_coeff(param_set)
    zi = SFP.gustiness_zi(param_set)
    g = SFP.grav(param_set)
    pv = SF.PointValueScheme()
    la = SF.LayerAverageScheme()

    @testset "Three-argument form" begin
        spec = SF.FlooredDeardorffGustinessSpec(FT(1))
        @test SF.gustiness_value(spec, param_set, FT(0)) == FT(1)
        @test SF.gustiness_value(spec, param_set, FT(-0.1)) == FT(1)
        B = FT(0.05)
        @test SF.gustiness_value(spec, param_set, B) ==
              max(FT(1), SF.gustiness_value(SF.DeardorffGustinessSpec(), param_set, B))
        @test SF.gustiness_value(SF.FlooredDeardorffGustinessSpec(FT(0)), param_set, B) ==
              SF.gustiness_value(SF.DeardorffGustinessSpec(), param_set, B)
        @test !SF.depends_on_ustar(spec)
    end

    @testset "Closed form at fixed ζ is a fixed point" begin
        spec = SF.FlooredDeardorffGustinessSpec(FT(0))
        for ΔT in (1, 10),
            ζ in (FT(-0.1), FT(-5)),
            z0m in (FT(0.01), FT(1)),
            sch in (pv, la)

            inp = make_inputs(FT, param_set, spec; ΔT, z0m)
            U = SF.free_convection_wind_speed(param_set, ζ, FT(0.3), inp, sch)
            @test U > 0
            # The buoyancy flux of the bulk relations at wind speed U, with the profile
            # integrals of the same scheme
            Δz_eff = SF.effective_height(inp)
            ϕ_m = SF.compute_physical_scale_coeff(
                param_set, Δz_eff, ζ, z0m, UF.MomentumTransport(), sch)
            ϕ_h = SF.compute_physical_scale_coeff(
                param_set, Δz_eff, ζ, z0m / 10, UF.HeatTransport(), sch)
            T_sfc = inp.T_sfc_guess
            ρ_sfc = SF.surface_density(
                param_set, inp.T_int, inp.ρ_int, T_sfc, inp.Δz, inp.q_tot_int,
                inp.q_liq_int, inp.q_ice_int, inp.q_vap_sfc_guess)
            θ_v_sfc, θ_v_int = SF.virtual_pottemps(
                param_set, inp, T_sfc, ρ_sfc, inp.q_vap_sfc_guess)
            ustar = U * ϕ_m
            θ_v_star = (θ_v_sfc - θ_v_int) * ϕ_h
            B = g / θ_v_int * ustar * θ_v_star
            @test β * cbrt(B * zi) ≈ U rtol = 1e-12
        end
        # The heat conductance evaluates the gustiness with its own scheme
        inp = make_inputs(FT, param_set, spec; U = 0)
        for sch in (pv, la)
            Ch = SF.heat_exchange_coefficient(
                param_set, FT(-1), FT(0.01), FT(0.001), SF.effective_height(inp), sch)
            @test SF.heat_conductance(
                param_set, FT(-1), FT(0.3), inp, FT(0.01), FT(0.001), sch) ≈
                  Ch * SF.gustiness_value(spec, param_set, FT(-1), FT(0.3), inp, sch)
        end
        # The schemes give different profile integrals, hence different closed forms
        inp = make_inputs(FT, param_set, spec; ΔT = 10)
        @test SF.free_convection_wind_speed(param_set, FT(-1), FT(0.3), inp, pv) !=
              SF.free_convection_wind_speed(param_set, FT(-1), FT(0.3), inp, la)
        # The scheme defaults to the point-value scheme
        @test SF.free_convection_wind_speed(param_set, FT(-1), FT(0.3), inp) ==
              SF.free_convection_wind_speed(param_set, FT(-1), FT(0.3), inp, pv)
        # No convective gustiness when the surface is not warmer than the air
        for ΔT in (0, -5)
            inp = make_inputs(FT, param_set, spec; ΔT)
            @test SF.free_convection_wind_speed(param_set, FT(0.1), FT(0.3), inp) == 0
            @test SF.gustiness_value(
                SF.FlooredDeardorffGustinessSpec(FT(0.5)), param_set, FT(0.1), FT(0.3), inp,
            ) == FT(0.5)
        end
        # The floor applies when it exceeds the convective part
        inp = make_inputs(FT, param_set, spec; ΔT = 0.01)
        @test SF.gustiness_value(
            SF.FlooredDeardorffGustinessSpec(FT(2)), param_set, FT(-0.01), FT(0.3), inp,
        ) == FT(2)
        # The friction velocity follows from ζ without an inner solve, with the
        # gustiness of the solver's scheme
        inp = make_inputs(FT, param_set, spec; U = 0)
        for sch in (pv, la)
            ustar, z0m_out, _ =
                @inferred SF.compute_ustar_and_roughness(param_set, FT(-1), inp, sch)
            @test ustar ≈
                  SF.gustiness_value(spec, param_set, FT(-1), ustar, inp, sch) *
                  SF.compute_physical_scale_coeff(
                param_set, SF.effective_height(inp), FT(-1), z0m_out,
                UF.MomentumTransport(), sch)
        end
    end

    @testset "Full solve" begin
        opts = SF.SolverOptions{FT}(maxiter = 40)
        for sch in (pv, la),
            u_min in (FT(0), FT(0.5)),
            U in (0, 0.1, 0.5, 2, 5),
            ΔT in (1, 10)

            spec = SF.FlooredDeardorffGustinessSpec(u_min)
            inp = make_inputs(FT, param_set, spec; U, ΔT)
            out = @inferred SF.surface_fluxes(param_set, inp, sch, opts)
            @test out.converged
            @test isfinite(out.shf) && 0 < out.shf < 2000
            @test out.ζ < 0
            @test FT(1e-3) < out.ustar < FT(4)
            ϕ_m = SF.compute_physical_scale_coeff(
                param_set, SF.effective_height(inp), out.ζ, FT(0.01),
                UF.MomentumTransport(), sch)
            # At the converged ζ, the closed form equals the Deardorff gustiness of the
            # buoyancy flux implied by ζ and u★, which the post-solve fluxes use. The
            # friction velocity the solve returns is therefore that of the post-solve
            # effective wind speed, to the convergence tolerance of ζ.
            B_ζ = SF.buoyancy_flux(param_set, out.ζ, out.ustar, inp)
            U_eff_ζ = max(FT(U), SF.gustiness_value(spec, param_set, B_ζ))
            @test out.ustar ≈ U_eff_ζ * ϕ_m rtol = 1e-6
            # The thermodynamic buoyancy flux of the returned heat fluxes agrees with
            # the similarity-theory value to the approximations relating θ_v★ to θ★
            # and q★
            ρ_sfc = SF.surface_density(
                param_set, inp.T_int, inp.ρ_int, inp.T_sfc_guess, inp.Δz, inp.q_tot_int,
                inp.q_liq_int, inp.q_ice_int, inp.q_vap_sfc_guess)
            B = SF.buoyancy_flux(
                param_set, out.shf, out.lhf, inp.T_sfc_guess, ρ_sfc, inp.q_vap_sfc_guess,
            )
            U_eff = max(FT(U), SF.gustiness_value(spec, param_set, B))
            @test out.ustar ≈ U_eff * ϕ_m rtol = 0.05
        end
        # Zero wind, zero floor: the free-convection limit gives finite, positive fluxes
        spec = SF.FlooredDeardorffGustinessSpec(FT(0))
        out0 =
            SF.surface_fluxes(param_set, make_inputs(FT, param_set, spec; U = 0), pv, opts)
        @test out0.converged && out0.shf > 10 && out0.ustar > 0.01
        # Stronger heating gives stronger convective exchange
        out1 = SF.surface_fluxes(
            param_set, make_inputs(FT, param_set, spec; U = 0, ΔT = 20), pv, opts)
        @test out1.shf > out0.shf
        # In stable conditions the result is that of the constant floor
        for (gust_a, gust_b) in (
            (SF.FlooredDeardorffGustinessSpec(FT(1)), SF.ConstantGustinessSpec(FT(1))),
            (SF.FlooredDeardorffGustinessSpec(FT(0)), SF.ConstantGustinessSpec(FT(0))),
        )
            a = SF.surface_fluxes(
                param_set, make_inputs(FT, param_set, gust_a; ΔT = -5), pv, opts)
            b = SF.surface_fluxes(
                param_set, make_inputs(FT, param_set, gust_b; ΔT = -5), pv, opts)
            @test a.shf ≈ b.shf rtol = 1e-6
            @test a.ustar ≈ b.ustar rtol = 1e-6
        end
        # More exchange than with the floor alone in unstable, calm conditions
        a = SF.surface_fluxes(
            param_set,
            make_inputs(FT, param_set, SF.FlooredDeardorffGustinessSpec(FT(1)); U = 0),
            pv, opts)
        b = SF.surface_fluxes(
            param_set,
            make_inputs(FT, param_set, SF.ConstantGustinessSpec(FT(1)); U = 0),
            pv, opts)
        @test a.shf > b.shf
    end

    @testset "Float32 and dual numbers" begin
        FT32 = Float32
        ps32 = SFP.SurfaceFluxesParameters(FT32, UF.BusingerParams)
        spec32 = SF.FlooredDeardorffGustinessSpec(FT32(1))
        inp32 = make_inputs(FT32, ps32, spec32; U = 0.5)
        out32 = @inferred SF.surface_fluxes(
            ps32, inp32, pv, SF.SolverOptions{FT32}(maxiter = 40))
        @test out32.shf isa FT32 && isfinite(out32.shf) && out32.shf > 0
        @test SF.free_convection_wind_speed(ps32, FT32(-1), FT32(0.3), inp32) isa FT32

        spec = SF.FlooredDeardorffGustinessSpec(FT(1))
        inp = make_inputs(FT, param_set, spec; U = 0.5)
        ζd = ForwardDiff.Dual(FT(-1), one(FT))
        Ud = @inferred SF.free_convection_wind_speed(param_set, ζd, FT(0.3), inp)
        @test Ud isa ForwardDiff.Dual
        @test ForwardDiff.value(Ud) ==
              SF.free_convection_wind_speed(param_set, FT(-1), FT(0.3), inp)
        ustar_d = (@inferred SF.compute_ustar_and_roughness(param_set, ζd, inp, pv))[1]
        @test ustar_d isa ForwardDiff.Dual
    end
end

end # module
