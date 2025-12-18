using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD
import ClimaParams

# Define a custom roughness model that uses LAI
struct LAIRoughnessParams{FT} <: SF.AbstractRoughnessParams
    base_z0::FT
end

# Define roughness methods for the custom model
function SF.momentum_roughness(
    spec::LAIRoughnessParams{FT},
    u★,
    sfc_param_set,
    roughness_inputs,
) where {FT}
    # Simple fake formula: z0 = base_z0 * LAI
    return spec.base_z0 * roughness_inputs.LAI
end

function SF.scalar_roughness(
    spec::LAIRoughnessParams{FT},
    u★,
    sfc_param_set,
    roughness_inputs,
) where {FT}
    return spec.base_z0 * roughness_inputs.LAI * FT(0.1)
end

function SF.momentum_and_scalar_roughness(
    spec::LAIRoughnessParams{FT},
    u★,
    sfc_param_set,
    roughness_inputs,
) where {FT}
    z0m = SF.momentum_roughness(spec, u★, sfc_param_set, roughness_inputs)
    z0s = SF.scalar_roughness(spec, u★, sfc_param_set, roughness_inputs)
    return (z0m, z0s)
end

@testset "Roughness Inputs Verification" begin
    FT = Float64
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
    thermo_params = SFP.thermodynamics_params(param_set)

    T_int = FT(300)
    p_int = FT(1e5)
    q_tot_int = FT(0.01)
    T_sfc_guess = FT(302) # Unstable
    q_vap_sfc_guess = FT(0.012)

    # Calculate density
    R_m = TD.gas_constant_air(thermo_params, q_tot_int, FT(0), FT(0))
    ρ_int = p_int / (R_m * T_int)

    # Custom configuration with our LAI model
    config = SF.SurfaceFluxConfig(
        LAIRoughnessParams(0.01),
        SF.gustiness_constant(1.0),
    )

    # Case 1: LAI = 1.0
    inputs1 = (LAI = 1.0,)
    result1 = SF.surface_fluxes(
        param_set,
        T_int,
        q_tot_int,
        ρ_int,
        T_sfc_guess,
        q_vap_sfc_guess,
        FT(0),
        FT(10),
        FT(0),
        nothing,
        nothing,
        inputs1, # roughness_inputs
        config,
    )

    # Case 2: LAI = 2.0 -> Higher roughness -> Higher Cd
    inputs2 = (LAI = 2.0,)
    result2 = SF.surface_fluxes(
        param_set,
        T_int,
        q_tot_int,
        ρ_int,
        T_sfc_guess,
        q_vap_sfc_guess,
        FT(0),
        FT(10),
        FT(0),
        nothing,
        nothing,
        inputs2, # roughness_inputs
        config,
    )

    @test result2.Cd > result1.Cd
    println("Cd (LAI=1): ", result1.Cd)
    println("Cd (LAI=2): ", result2.Cd)
end
