# Maximum relative mismatch tolerated between the point-value and layer-average
# velocity-scale computations. With the recovery profiles below, the schemes
# usually agree to within ≲ 15%, so this limit still flags meaningful
# regressions while allowing the expected discretization differences.
const RECOVERY_SCHEME_RTOL = FloatType(0.15)

# Helper to compute the physical scale coefficient with allocation checks.
# This is used in the "SurfaceFluxes - Finite-Difference vs Finite-Volume" testset.
function compute_physical_scale_coeff_with_checks(
    param_set,
    sc,
    L_MO,
    scheme,
    z0m,
    z0b,
)
    @test_allocs_and_ts SF.compute_physical_scale_coeff(
        param_set,
        sc,
        L_MO,
        UF.MomentumTransport(),
        scheme,
        z0m,
        z0b,
    )
end

@testset "SurfaceFluxes - Finite-Difference vs Finite-Volume" begin
    # Verify that different scale computation schemes (point value vs layer
    # average) produce consistent velocity scales for a set of representative 
    # test cases and that the various surface-condition containers agree for the
    # same discretization scheme.

    param_set = SFP.SurfaceFluxesParameters(FloatType, BusingerParams)
    thermo_params = param_set.thermo_params

    ρ_sfc = FloatType(1.15)
    ρ_in = FloatType(1.13)
    qt_sfc = FloatType(0.01)
    qt_in = FloatType(0.009)

    # Discretization altitude z [m]
    z = ArrayType(
        FloatType[
            29.432779269303,
            30.0497139076724,
            31.6880000418153,
            34.1873479240475,
        ],
    )

    # Surface temperature [K]
    T_sfc = ArrayType(
        FloatType[
            277.38348,
            276.8474,
            280.4166,
            296.93747,
        ],
    )

    # Temperature at interior level at height z [K]
    T_in = ArrayType(
        FloatType[
            271.62735,
            272.46533,
            277.3738,
            291.56274,
        ],
    )

    # Roughness lengths for momentum [m]
    z0 = ArrayType(
        FloatType[
            5.86144925739178e-05,
            0.0001,
            0.000641655193293549,
            3.23383768877187e-05,
        ],
    )

    # Wind speed at z [m s⁻¹]
    speed = ArrayType(
        FloatType[
            2.9693638452068,
            2.43308757772094,
            5.69418282305367,
            9.5608693754561,
        ],
    )

    # Friction velocity [m s⁻¹]
    u_star = ArrayType(
        FloatType[
            0.109462510724615,
            0.0932942802513508,
            0.223232887323184,
            0.290918439028557,
        ],
    )

    # Target Monin–Obukhov lengths [m]; the set spans unstable through very
    # stable surface layers and keeps the point-value vs layer-average schemes
    # within ≲ 15% of each other.
    L_MO_targets = ArrayType(
        FloatType[
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

            # TODO: Remove thermodynamic states from here and use the signature 
            # Create thermodynamic states from temperature, density, and specific humidity
            ts_sfc = TD.PhaseEquil_ρTq(thermo_params, ρ_sfc, T_sfc[ii], qt_sfc)
            ts_in = TD.PhaseEquil_ρTq(thermo_params, ρ_in, T_in[ii], qt_in)

            state_sfc =
                SF.StateValues(FloatType(0), (FloatType(0), FloatType(0)), ts_sfc)
            state_in =
                SF.StateValues(z[ii], (FloatType(speed[ii]), FloatType(0)), ts_in)

            # State containers for different computation modes
            z0m = z0[ii]
            z0b = FloatType(0.001)

            state_containers = (
                SF.Fluxes(
                    state_in,
                    state_sfc,
                    FloatType(0),
                    FloatType(0),
                    z0m,
                    z0b,
                ),
                SF.FluxesAndFrictionVelocity(
                    state_in,
                    state_sfc,
                    FloatType(0),
                    FloatType(0),
                    u_star[ii],
                    z0m,
                    z0b,
                ),
                SF.ValuesOnly(state_in, state_sfc, z0m, z0b),
            )

            point_scales = FloatType[]
            layer_scales = FloatType[]

            for sc in state_containers
                # Point-value scheme
                u_scale_fd = compute_physical_scale_coeff_with_checks(
                    param_set,
                    sc,
                    L_MO,
                    SF.PointValueScheme(),
                    z0m,
                    z0b,
                )

                # Layer-averaged scheme
                u_scale_fv = compute_physical_scale_coeff_with_checks(
                    param_set,
                    sc,
                    L_MO,
                    SF.LayerAverageScheme(),
                    z0m,
                    z0b,
                )

                @test isfinite(u_scale_fd) && u_scale_fd > FloatType(0)
                @test isfinite(u_scale_fv) && u_scale_fv > FloatType(0)

                rel_diff = abs(u_scale_fd - u_scale_fv) / u_scale_fd
                @test rel_diff <= RECOVERY_SCHEME_RTOL

                push!(point_scales, u_scale_fd)
                push!(layer_scales, u_scale_fv)
            end

            # Different surface-condition containers should give identical
            # velocity scales for a fixed discretization scheme.
            for scales in (point_scales, layer_scales)
                reference = first(scales)
                for candidate in scales
                    @test isapprox(candidate, reference; rtol = FloatType(5e-6))
                end
            end
        end
    end
end
