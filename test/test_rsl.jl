using Plots
using LaTeXStrings

import SurfaceFluxes 
const SF = SurfaceFluxes
SurfaceFluxes.error_on_non_convergence() = true

import SurfaceFluxes.UniversalFunctions as UF
using StaticArrays
using Thermodynamics
Thermodynamics.print_warning() = false
using BenchmarkTools
using RootSolvers
const RS = RootSolvers
include(joinpath(pkgdir(SurfaceFluxes), "parameters", "create_parameters.jl"))
const SFP = SF.Parameters

const FT = Float64;

### Generate parameter lists
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
param_set = create_parameters(toml_dict, UF.Gryanik())
thermo_params = SFP.thermodynamics_params(param_set)
uft = SFP.universal_func_type(param_set)
### 


# TODO: Implement quantitative tests for RSL/Non-RSL profiles