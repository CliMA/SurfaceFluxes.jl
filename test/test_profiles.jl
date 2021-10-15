using CLIMAParameters: AbstractEarthParameterSet
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()
FT = Float64

args = (param_set, [0.1, 0.1, 260.0], [0.0, 265.0], [0.0, 265.0], 0.1, 265.0, 25.0)
(EarthParameterSet(), [0.1, 0.1, 260.0], [8.0, 265.0], [0.0, 265.0], 0.1, 265.0, 25.0)
LMO_init = FT(0.1)
u_star_init = FT(0.1)
th_star_init = FT(260)
using SurfaceFluxes
result = surface_conditions(args..., SurfaceFluxes.FVScheme())
