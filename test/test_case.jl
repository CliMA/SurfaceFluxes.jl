import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD
import Thermodynamics.Parameters as TDP
import ClimaParams as CP

FT = Float32

param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)
thermo_params = SFP.thermodynamics_params(param_set)

ρ_air = FT(1.2)
T_air = FT(282)
T_sfc = FT(320)

roughness = SF.ConstantRoughnessParams(FT(1e-2), FT(1e-2))
gustiness_spec = SF.ConstantGustinessSpec(FT(1))
config = SF.SurfaceFluxConfig(roughness, gustiness_spec)
flux_specs = SF.FluxSpecs{FT}()

output = SF.surface_fluxes(
    param_set,
    T_air,
    FT(0),
    FT(0),
    FT(0),
    ρ_air,
    T_sfc,
    FT(0),
    FT(0),
    32,
    0,
    (FT(3), FT(0)),
    (FT(0), FT(0)),
    nothing,
    config, UF.PointValueScheme(), SF.SolverOptions{FT}(1e-2, 15, true)
)
