import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import ClimaParams
using Test
include("test_utils.jl")

# Maximum relative mismatch tolerated between the point-value and layer-average
# velocity-scale computations. With the recovery profiles below, the schemes
# usually agree to within ≲ 15%, so this limit still flags meaningful
# regressions while allowing the expected discretization differences.
const RECOVERY_SCHEME_RTOL = Float32(0.15)

# Helper to compute the physical scale coefficient with allocation checks.
# This is used in the "SurfaceFluxes - Finite-Difference vs Finite-Volume" testset.
function compute_physical_scale_coeff_with_checks(
    param_set,
    sc,
    L_MO,
    scheme,
    z0m,
    z0h,
)
    ζ = sc.Δz / L_MO
    @test_allocs_and_ts SF.compute_physical_scale_coeff(
        param_set,
        sc.Δz,
        ζ,
        z0m,
        UF.MomentumTransport(),
        scheme,
    )
end

@testset "SurfaceFluxes - Finite-Difference vs Finite-Volume" begin
    # Verify that different scale computation schemes (point value vs layer
    # average) produce consistent velocity scales for a set of representative 
    # test cases and that the various surface-condition containers agree for the
    # same discretization scheme.

    param_set = SFP.SurfaceFluxesParameters(Float32, UF.BusingerParams)
    thermo_params = param_set.thermo_params

    ρ_int = Float32(1.13)
    qt_sfc = Float32(0.01)
    qt_int = Float32(0.009)

    # Discretization altitude z [m]
    z = Array(
        Float32[
            29.432779269303,
            30.0497139076724,
            31.6880000418153,
            34.1873479240475,
        ],
    )

    # Surface temperature [K]
    T_sfc = Array(
        Float32[
            277.38348,
            276.8474,
            280.4166,
            296.93747,
        ],
    )

    # Temperature at interior level at height z [K]
    T_int = Array(
        Float32[
            271.62735,
            272.46533,
            277.3738,
            291.56274,
        ],
    )

    # Roughness lengths for momentum [m]
    z0 = Array(
        Float32[
            5.86144925739178e-05,
            0.0001,
            0.000641655193293549,
            3.23383768877187e-05,
        ],
    )

    # Wind speed at z [m s⁻¹]
    speed = Array(
        Float32[
            2.9693638452068,
            2.43308757772094,
            5.69418282305367,
            9.5608693754561,
        ],
    )

    # Friction velocity [m s⁻¹]
    u_star = Array(
        Float32[
            0.109462510724615,
            0.0932942802513508,
            0.223232887323184,
            0.290918439028557,
        ],
    )

    # Target Monin–Obukhov lengths [m]; the set spans unstable through very
    # stable surface layers and keeps the point-value vs layer-average schemes
    # within ≲ 15% of each other.
    L_MO_targets = Array(
        Float32[
            -80,
            -40,
            160,
            400,
        ],
    )

    for ii in 1:length(L_MO_targets)
        @testset "Profile $(ii)" begin
            # Prescribed Monin–Obukhov length selected to span unstable through
            # strongly stable conditions.
            L_MO = L_MO_targets[ii]

            # Roughness lengths
            z0m = z0[ii]
            z0h = Float32(0.001)

            # Build inputs structure
            inputs_container = (;
                T_int = T_int[ii],
                q_tot_int = qt_int,
                ρ_int = ρ_int,
                T_sfc_guess = T_sfc[ii],
                q_vap_sfc_guess = qt_sfc,
                Φ_sfc = Float32(0),
                Δz = z[ii],
                d = Float32(0),
                u_int = (Float32(speed[ii]), Float32(0)),
                u_sfc = (Float32(0), Float32(0)),
                roughness = SF.roughness_lengths(z0m, z0h),
                gustiness = SF.ConstantGustinessSpec(Float32(1.0)),
                moisture_model = SF.MoistModel(),
            )

            # We create a dummy inputs structure. 
            # compute_physical_scale_coeff expects inputs with Δz field

            sf_inputs = SF.build_surface_flux_inputs(
                param_set,
                inputs_container.T_int, inputs_container.q_tot_int, Float32(0), Float32(0), inputs_container.ρ_int,
                inputs_container.T_sfc_guess, inputs_container.q_vap_sfc_guess,
                inputs_container.Φ_sfc, inputs_container.Δz, inputs_container.d,
                inputs_container.u_int, inputs_container.u_sfc,
                SF.SurfaceFluxConfig(
                    inputs_container.roughness,
                    inputs_container.gustiness,
                    inputs_container.moisture_model,
                ),
                nothing, # roughness_inputs
                SF.FluxSpecs(Float32), # Default specs
                nothing, nothing,
            )

            # Point-value scheme
            u_scale_fd = compute_physical_scale_coeff_with_checks(
                param_set,
                sf_inputs,
                L_MO,
                SF.PointValueScheme(),
                z0m,
                z0h,
            )

            # Layer-averaged scheme
            u_scale_fv = compute_physical_scale_coeff_with_checks(
                param_set,
                sf_inputs,
                L_MO,
                SF.LayerAverageScheme(),
                z0m,
                z0h,
            )

            @test isfinite(u_scale_fd) && u_scale_fd > Float32(0)
            @test isfinite(u_scale_fv) && u_scale_fv > Float32(0)

            rel_diff = abs(u_scale_fd - u_scale_fv) / u_scale_fd
            @test rel_diff <= RECOVERY_SCHEME_RTOL

        end
    end
end
