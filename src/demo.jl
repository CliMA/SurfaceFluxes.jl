import SurfaceFluxes
SurfaceFluxes.error_on_non_convergence() = true
const SF = SurfaceFluxes
import SurfaceFluxes.UniversalFunctions as UF
using StaticArrays
using Thermodynamics
using Thermodynamics.TemperatureProfiles
using Thermodynamics.TestedProfiles
Thermodynamics.print_warning() = false
using BenchmarkTools

using RootSolvers
const RS = RootSolvers

include(joinpath(pkgdir(SurfaceFluxes), "parameters", "create_parameters.jl"))
const SFP = SF.Parameters
import Thermodynamics.TestedProfiles: input_config, PhaseEquilProfiles


FT = Float32
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
param_set = create_parameters(toml_dict, UF.Businger())
thermo_params = SFP.thermodynamics_params(param_set)
z_int = FT(250)
z_sfc = FT(0)
ts_int = Thermodynamics.PhaseEquil{FT}(1.1749369, 97086.64, 10584.462, 0.0, 287.91174) 
ts_sfc = Thermodynamics.PhaseEquil{FT}(1.1575673, 95651.26, 32775.72, 0.009824766, 286.2016)
z0m, z0b = FT(1.0f-5), FT(1.0f-5)
sc = SF.ValuesOnly{FT}(;
                       state_in = SF.InteriorValues(z_int, (FT(0.5), FT(0.5)), ts_int),
                       state_sfc = SF.SurfaceValues(z_sfc, (FT(0), FT(0)), ts_sfc),
                       z0m = z0m,
                       z0b = z0b,
                       )
gryanik_analytical_fd = SF.surface_conditions(param_set, sc, SF.FDScheme(); noniterative_stable_sol = true)
classical_iterative_fd = SF.surface_conditions(param_set, sc, SF.FDScheme(); noniterative_stable_sol = false)
gryanik_analytical_fv = SF.surface_conditions(param_set, sc, SF.FVScheme(); noniterative_stable_sol = true)
classical_iterative_fv = SF.surface_conditions(param_set, sc, SF.FVScheme(); noniterative_stable_sol = false)
@info "Analytical Solution FD" gryanik_analytical_fd;
@info "Standard Solution FD" classical_iterative_fd;
@info "Analytical Solution FV" gryanik_analytical_fv;
@info "Standard Solution FV" classical_iterative_fv;
