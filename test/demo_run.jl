import Random
Random.seed!(1234)
import Thermodynamics as TD
using SurfaceFluxes
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions.BusingerParams
import ClimaParams as CP
import LinearAlgebra: norm
using SurfaceFluxes: SimilarityScaleVars

ArrayType = Array
FT = Float32

@info "Creating Input Struct"
u★ = FT(0.1)
DSEᵥ★ = FT(10)
q★ = FT(0.0001)
L★ = FT(-1)
𝓁u = FT(0.001)
𝓁θ = FT(0.0001)
@info "Success @ Creating Input Struct"


@info "Setting up inputs"
param_set = SFP.SurfaceFluxesParameters(FT, BusingerParams)
thermo_params = param_set.thermo_params
uft = UF.universal_func_type(typeof(param_set.ufp))
ρ_sfc = FT(1.15)
ρ_in = FT(1.13)
qt_sfc = FT(0.01)
qt_in = FT(0.009)
z = FT(29.432779269303)
## Virtual Potential Temperature at height z
θ = FT(268.559120403867)
## Surface Pottemp
θ_sfc = FT(273.42369841804)
## Roughness lengths
z0 = FT(0.0001)
speed = FT(2.9693638452068)
## Scale velocity and moisture
u★ = FT(0.109462510724615)
b★ = FT(-0.00690834676781433)
𝜅 = SFP.von_karman_const(param_set)
L_MO = u★^2/𝜅/b★

ts_sfc = TD.PhaseEquil_ρθq(thermo_params, ρ_sfc, θ_sfc, qt_sfc)
ts_in = TD.PhaseEquil_ρθq(thermo_params, ρ_in, θ, qt_in)
state_sfc =
    SF.StateValues(FT(0), (FT(0), FT(0)), ts_sfc)
state_in =
    SF.StateValues(z, (FT(speed), FT(0)), ts_in)

z0m = z0
z0b = FT(0.00011)

sc = SF.ValuesOnly(state_in, state_sfc, z0m, z0b, roughness_model = SF.CharnockRoughness())
#sc = SF.ValuesOnly(state_in, state_sfc, z0m, z0b, roughness_model=SF.ScalarRoughness())
#SF.obukhov_iteration(X★, sc,uft, scheme, param_set)
