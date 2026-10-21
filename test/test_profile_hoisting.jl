# The stability solve evaluates the dimensionless profiles once per residual
# evaluation and stores the coefficients `ϕ_m` and `ϕ_h` in the inputs
# (`with_profile_coefficients`). The tests check that the stored coefficients
# equal those of `compute_physical_scale_coeff`, that the readers return them,
# and that the heat conductance and the fluxes of the solve match those of
# inputs without stored coefficients.

module TestProfileHoisting

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import ClimaParams as CP
using Main: @test_allocs_and_ts

# Allocation and type-stability checks need the arguments passed in, since
# locals of a testset body captured by the check's closure are boxed
with_coefficients_checked(inputs, ps, scheme, ζ) =
    @test_allocs_and_ts SF.with_profile_coefficients(
        inputs,
        ps,
        scheme,
        ζ,
        inputs.z0m,
        inputs.z0h,
    )

@testset "Profile coefficients" begin
    for FT in (Float32, Float64), scheme in (SF.PointValueScheme(), SF.LayerAverageScheme())
        param_set = SFP.SurfaceFluxesParameters(FT, UF.GryanikParams)
        config = SF.SurfaceFluxConfig(
            SF.ConstantRoughnessParams{FT}(FT(0.1), FT(0.01)),
            SF.FlooredDeardorffGustinessSpec(FT(0.5)),
            SF.MoistModel(),
            SF.NoRoughnessSubLayer(),
            SF.MaxHeatFluxStabilityCap(),
            SF.ReferenceAboveSurface(),
        )
        bare = SF.build_surface_flux_inputs(
            FT(290), FT(0.005), FT(0), FT(0), FT(1.2), FT(295), FT(0.012), FT(0),
            FT(20),
            FT(5), (FT(1), FT(0)), (FT(0), FT(0)), config, nothing, SF.FluxSpecs{FT}(),
            nothing, nothing,
        )
        inputs = SF.with_stability_cap(
            SF.with_roughness_lengths(bare, param_set),
            param_set,
            scheme,
        )
        @test inputs.ϕ_m === nothing && inputs.ϕ_h === nothing
        for ζ in (FT(-2), FT(-0.1), FT(0), FT(0.3), FT(5))
            stored = with_coefficients_checked(inputs, param_set, scheme, ζ)
            Δz_eff = SF.effective_height(param_set, inputs)
            ζ_c = SF.capped_stability(inputs, ζ)
            ϕ_m = SF.compute_physical_scale_coeff(
                param_set, Δz_eff, ζ_c, inputs.z0m, UF.MomentumTransport(), scheme,
            )
            ϕ_h = SF.compute_physical_scale_coeff(
                param_set, Δz_eff, ζ_c, inputs.z0h, UF.HeatTransport(), scheme,
            )
            @test stored.ϕ_m === ϕ_m && stored.ϕ_h === ϕ_h
            # Readers return the stored coefficients when present
            @test (@inferred SF.physical_scale_coeff(
                param_set, stored, scheme, ζ, inputs.z0m, UF.MomentumTransport(),
            )) === ϕ_m
            @test SF.physical_scale_coeff(
                param_set, inputs, scheme, ζ, inputs.z0h, UF.HeatTransport(),
            ) === ϕ_h
            # Heat conductance from the solver's friction velocity equals Ch ΔU
            ΔU = SF.windspeed(param_set, ζ, zero(FT), inputs, scheme)
            u★ = ΔU * ϕ_m
            g_h = SF.heat_conductance(
                param_set,
                ζ,
                u★,
                stored,
                inputs.z0m,
                inputs.z0h,
                scheme,
            )
            @test g_h ≈ SF.heat_conductance(
                param_set,
                ζ,
                u★,
                inputs,
                inputs.z0m,
                inputs.z0h,
                scheme,
            )
            @test g_h ≈ ϕ_m * ϕ_h * ΔU
            # The bulk Richardson number of the profiles
            @test SF.bulk_richardson_number(
                param_set,
                stored,
                scheme,
                ζ,
                inputs.z0m,
                inputs.z0h,
            ) ≈
                  SF.bulk_richardson_number(
                SFP.uf_params(param_set), inputs.rsl_model, Δz_eff, ζ, inputs.z0m,
                inputs.z0h, scheme, inputs.ζ_cap,
            )
        end
        # The interior virtual potential temperature stored once per solve
        hoisted = SF.with_interior_virtual_pottemp(inputs, param_set)
        @test hoisted.θ_v_int === SF.interior_virtual_pottemp(param_set, inputs)
        @test SF.virtual_pottemps(param_set, hoisted, FT(295), FT(1.2), FT(0.01)) ===
              SF.virtual_pottemps(param_set, inputs, FT(295), FT(1.2), FT(0.01))

        # Inside a solve, callbacks see the stored coefficients for roughness
        # lengths independent of u★, and θ_v_int for both roughness models
        seen = Ref(true)
        function record_T(ζ, ps, thermo_params, inputs, scheme, u★, z0m, z0h)
            seen[] &= (inputs.ϕ_m !== nothing) & (inputs.ϕ_h !== nothing)
            seen[] &= inputs.θ_v_int !== nothing
            return inputs.T_sfc_guess
        end
        function record_T_coare(ζ, ps, thermo_params, inputs, scheme, u★, z0m, z0h)
            seen[] &= inputs.ϕ_m === nothing
            seen[] &= inputs.θ_v_int !== nothing
            return inputs.T_sfc_guess
        end
        for (roughness, record) in (
            (SF.ConstantRoughnessParams{FT}(FT(0.1), FT(0.01)), record_T),
            (SF.COARE3RoughnessParams{FT}(), record_T_coare),
        )
            seen[] = true
            config_r = SF.SurfaceFluxConfig(
                roughness, SF.FlooredDeardorffGustinessSpec(FT(0.5)), SF.MoistModel(),
                SF.NoRoughnessSubLayer(), SF.MaxHeatFluxStabilityCap(),
                SF.ReferenceAboveSurface(),
            )
            sc = SF.surface_fluxes(
                param_set, FT(290), FT(0.005), FT(0), FT(0), FT(1.2), FT(295),
                FT(0.012),
                FT(0), FT(20), FT(5), (FT(1), FT(0)), (FT(0), FT(0)), nothing, config_r,
                scheme, nothing, nothing, record, nothing,
            )
            @test sc.converged && seen[]
        end
    end
end

end # module
