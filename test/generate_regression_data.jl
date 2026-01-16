import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import ClimaParams as CP
using Printf

function generate_cases()
    FT = Float32
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)

    # --- GRID DEFINITION ---
    u_vals = [FT(0.1), FT(5.0), FT(20.0)]
    dTheta_vals = [FT(5.0), FT(0.1), FT(0.0), FT(-0.1), FT(-5.0)] # Stable to Unstable
    dq_vals = [FT(-0.002), FT(0.0), FT(0.002)]

    roughness_configs = [
        (name = "Fixed_0.0001", model = SF.ConstantRoughnessParams(FT(0.0001), FT(0.0001))),
        (name = "Fixed_0.1", model = SF.ConstantRoughnessParams(FT(0.1), FT(0.1))),
        (name = "Fixed_Default", model = SF.ConstantRoughnessParams(FT(1e-3), FT(1e-3))),
        (name = "COARE3", model = SF.COARE3RoughnessParams{FT}()),
    ]

    T_sfc = FT(300.0)
    q_sfc = FT(0.015)

    Δz = FT(10.0)
    ρ_int = FT(1.2)

    maxiter = 100
    tol = FT(1e-10)
    opts = SF.SolverOptions{FT}(maxiter = maxiter, tol = tol, rtol = tol)

    cases = []

    count = 1
    for u in u_vals
        for dTheta in dTheta_vals
            for dq in dq_vals
                for r_config in roughness_configs

                    T_air = T_sfc + dTheta
                    q_air = q_sfc + dq

                    u_int_val = (u, FT(0))
                    u_sfc_val = (FT(0), FT(0))

                    roughness_inputs = nothing
                    gustiness = SF.ConstantGustinessSpec(FT(1.0))

                    config =
                        SF.SurfaceFluxConfig(r_config.model, gustiness, SF.MoistModel())

                    inputs = SF.build_surface_flux_inputs(
                        T_air, q_air, FT(0), FT(0), ρ_int,
                        T_sfc, q_sfc,
                        FT(0), Δz, FT(0),
                        u_int_val, u_sfc_val,
                        config,
                        roughness_inputs,
                        SF.FluxSpecs{FT}(),
                        nothing, nothing,
                    )

                    result =
                        SF.surface_fluxes(param_set, inputs, SF.PointValueScheme(), opts)

                    case_val = (;
                        name = "Case $count: u=$(u), dT=$(dTheta), dq=$(dq), $(r_config.name)",
                        u_int = u_int_val,
                        T_int = T_air,
                        q_tot_int = q_air,
                        ρ_int = ρ_int,
                        T_sfc = T_sfc,
                        q_sfc = q_sfc,
                        config_name = r_config.name,
                        expected = (;
                            shf = result.shf,
                            lhf = result.lhf,
                            ustar = result.ustar,
                            Cd = result.Cd,
                            g_h = result.g_h,
                            evaporation = result.evaporation,
                        ),
                    )

                    push!(cases, case_val)
                    count += 1
                end
            end
        end
    end

    return cases
end

function print_cases(cases)
    println("function case_definitions(::Type{FT}) where {FT}")
    println("    return [")

    for c in cases
        r_code = ""
        # Match names to construct valid code
        if c.config_name == "Fixed_0.0001"
            r_code = "SF.ConstantRoughnessParams(FT(0.0001), FT(0.0001))"
        elseif c.config_name == "Fixed_0.1"
            r_code = "SF.ConstantRoughnessParams(FT(0.1), FT(0.1))"
        elseif c.config_name == "Fixed_Default"
            r_code = "SF.ConstantRoughnessParams(FT(1e-3), FT(1e-3))"
        elseif c.config_name == "COARE3"
            r_code = "SF.COARE3RoughnessParams{FT}(z0m_default=FT(1e-3))"
        end

        # Use repr to ensure float precision in output
        println("    (")
        println("        name = \"$(c.name)\",")
        println("        u_int = (FT($(c.u_int[1])), FT($(c.u_int[2]))),")
        println("        u_sfc = (FT(0), FT(0)),")
        println("        T_int = FT($(c.T_int)),")
        println("        q_tot_int = FT($(c.q_tot_int)),")
        println("        ρ_int = FT($(c.ρ_int)),")
        println("        T_sfc = FT($(c.T_sfc)),")
        println("        q_sfc = FT($(c.q_sfc)),")
        println("        Δz = FT(10.0),")
        println("        roughness_config = $r_code,")
        println("        expected = (;")
        println("            shf = FT($(c.expected.shf)),")
        println("            lhf = FT($(c.expected.lhf)),")
        println("            ustar = FT($(c.expected.ustar)),")
        println("            Cd = FT($(c.expected.Cd)),")
        println("            g_h = FT($(c.expected.g_h)),")
        println("            evaporation = FT($(c.expected.evaporation)),")
        println("        ),")
        println("    ),")
    end
    println("    ]")
    println("end")
end

cases = generate_cases()
print_cases(cases)
