# A user-defined roughness model receives the `roughness_inputs` of a solve: here, a
# roughness length proportional to a plant area index `PAI`, and a solve whose drag
# coefficient increases with it. The model does not declare `depends_on_ustar`, so the
# solver treats it as dependent on the friction velocity.

module TestRoughnessInputs

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD
import ClimaParams

struct PAIRoughnessParams{FT} <: SF.AbstractRoughnessParams
    base_z0::FT
end

# Define roughness methods for the custom model
function SF.momentum_roughness(
    spec::PAIRoughnessParams{FT},
    u★,
    sfc_param_set,
    roughness_inputs,
) where {FT}
    return spec.base_z0 * roughness_inputs.PAI
end

function SF.scalar_roughness(
    spec::PAIRoughnessParams{FT},
    u★,
    sfc_param_set,
    roughness_inputs,
) where {FT}
    return spec.base_z0 * roughness_inputs.PAI * FT(0.1)
end

function SF.momentum_and_scalar_roughness(
    spec::PAIRoughnessParams{FT},
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

    # Custom configuration with the PAI model
    config = SF.SurfaceFluxConfig(
        PAIRoughnessParams(0.01),
        SF.ConstantGustinessSpec(1.0),
    )

    u_int, u_sfc = (FT(3), FT(0)), (FT(0), FT(0))
    inputs1 = (PAI = FT(1),)
    result1 = SF.surface_fluxes(
        param_set,
        T_int,
        q_tot_int,
        FT(0),
        FT(0),
        ρ_int,
        T_sfc_guess,
        q_vap_sfc_guess,
        FT(0),
        FT(10),
        FT(0),
        u_int,
        u_sfc,
        inputs1, # roughness_inputs
        config,
    )

    inputs2 = (PAI = FT(2),)
    result2 = SF.surface_fluxes(
        param_set,
        T_int,
        q_tot_int,
        FT(0),
        FT(0),
        ρ_int,
        T_sfc_guess,
        q_vap_sfc_guess,
        FT(0),
        FT(10),
        FT(0),
        u_int,
        u_sfc,
        inputs2, # roughness_inputs
        config,
    )

    @test result1.converged && result2.converged
    # A rougher surface has a larger drag coefficient
    @test result2.Cd > result1.Cd
end

end # module
