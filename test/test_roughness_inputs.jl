using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD
import ClimaParams

# Define a custom roughness model that uses LAI
struct LAIRoughnessModel{FT} <: SF.AbstractRoughnessModel{FT}
    base_z0::FT
end

# Define roughness methods for the custom model
# Note: we need to import these to extend them if they are not exported, 
# or just define them if we are in the test module and SF exports them.
# SF exports them? Let's check. No, they are not exported.
# So we need to extend SF.momentum_roughness etc.

function SF.momentum_roughness(model::LAIRoughnessModel{FT}, u★, sfc_param_set, ctx, roughness_inputs) where {FT}
    # Simple fake formula: z0 = base_z0 * LAI
    return model.base_z0 * roughness_inputs.LAI
end

function SF.scalar_roughness(model::LAIRoughnessModel{FT}, u★, sfc_param_set, ctx, roughness_inputs) where {FT}
    return model.base_z0 * roughness_inputs.LAI * FT(0.1)
end

@testset "Roughness Inputs Verification" begin
    FT = Float64
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
    thermo_params = SFP.thermodynamics_params(param_set)

    Tin = FT(300)
    qin = FT(0.01)
    Ts_guess = FT(302) # Unstable
    qs_guess = FT(0.012)
    
    # Calculate density
    q_pt = TD.PhasePartition_equil_given_p(
        thermo_params,
        Tin,
        FT(100000),
        qin,
        TD.PhaseEquil,
    )
    ρin = TD.air_density(thermo_params, Tin, FT(100000), q_pt)

    # Custom configuration with our LAI model
    # We need to bypass the config struct if it enforces types, or make a custom spec.
    # SurfaceFluxConfig requires AbstractRoughnessSpec.
    
    struct LAIRoughnessSpec{FT} <: SF.AbstractRoughnessSpec
        base_z0::FT
    end
    
    # We need to extend instantiate_roughness
    function SF.instantiate_roughness(param_set::SF.APS{FT}, spec::LAIRoughnessSpec) where {FT}
        return LAIRoughnessModel{FT}(convert(FT, spec.base_z0))
    end

    config = SF.SurfaceFluxConfig(
        LAIRoughnessSpec(0.01),
        SF.gustiness_constant(1.0)
    )

    # Case 1: LAI = 1.0
    inputs1 = (LAI = 1.0,)
    result1 = SF.surface_fluxes(
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
        inputs1, # roughness_inputs
        config,
    )

    # Case 2: LAI = 2.0 -> Higher roughness -> Higher Cd
    inputs2 = (LAI = 2.0,)
    result2 = SF.surface_fluxes(
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
        inputs2, # roughness_inputs
        config,
    )

    @test result2.Cd > result1.Cd
    println("Cd (LAI=1): ", result1.Cd)
    println("Cd (LAI=2): ", result2.Cd)
end
