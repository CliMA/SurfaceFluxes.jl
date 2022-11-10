using NCDatasets

using SurfaceFluxes
const SF = SurfaceFluxes
SurfaceFluxes.error_on_non_convergence() = true
import SurfaceFluxes.UniversalFunctions as UF
using StaticArrays
using Thermodynamics
using Thermodynamics.TemperatureProfiles
using Thermodynamics.TestedProfiles
Thermodynamics.print_warning() = false

using RootSolvers
const RS = RootSolvers

include(joinpath(pkgdir(SurfaceFluxes), "parameters", "create_parameters.jl"))
const SFP = SF.Parameters
import Thermodynamics.TestedProfiles: input_config, PhaseEquilProfiles

# Generate Parameter Lists
FT = Float32
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
surface_params = create_parameters(toml_dict, UF.Businger())
thermo_params = SFP.thermodynamics_params(surface_params)

# 950hPa level
extract_var(ds, n; 
            ilat= 1, 
            ilon= 1, 
            it= 1) = ds[n][ilat,ilon,it] # * ds[n].attrib["scale_factor"] + ds[n].attrib["add_offset"] # NCDatasets deals with this

data_1 = NCDataset("test_global.nc","r") do ds
    T = extract_var(ds, "t") # temperature
    u = extract_var(ds, "u") # eastward velo
    v = extract_var(ds, "v") # northward velo
    z = extract_var(ds, "z")/ 9.81 # altitude from "z" geopotential
    q = extract_var(ds, "q") # total moisture (specific humidity)
    p = 95000.0 # source data pressure level (950hPa)
    ts = TD.PhaseEquil_pTq(thermo_params, FT(p), FT(T), FT(q))
    ρ = air_density(thermo_params, ts)
    (; T = T, u = FT(u), v = FT(v), z = FT(z), q = q, ρ = ρ, 
    ts=ts)
end

# 2-m / surface level
data_sfc = NCDataset("test_global_sfc.nc","r") do ds
    T = extract_var(ds, "t2m") # temperature
    u = 0.0 # eastward velo (may be non zero for ocean)
    v = 0.0 # northward velo
    z = 0.0 # altitude from "z" geopotential

    e = saturation_vapor_pressure(thermo_params, T, TD.Liquid())
    # tetens_eq(𝒯) = 0.61078 * exp(17.27 * 𝒯  / (𝒯 + 237.3)) 
    # e = tetens_eq(T)
    # TODO: Update with ClimaParameters values 
    # (EarthParameters <: AbstractParameterSet)
    R = 287.1
    R_m = 461.5
    p = extract_var(ds, "sp")

    q_sat = R / R_m * e * (p - ((1 - R / R_m) * e)) # saturated total moisture (specific humidity)

    ts = TD.PhaseEquil_pTq(thermo_params, FT(p), FT(T), FT(q_sat))
    ρ = air_density(thermo_params, ts)
    # T_dew = extract_var(ds, "d2m"][:] # dewpoint temperature
    
    (; shf=extract_var(ds, "sshf"), lhf=extract_var(ds, "slhf"), evap = extract_var(ds, "e"),  
    T = T, u = FT(u), v = FT(v), z = FT(z), q = q_sat, ρ = ρ, 
    ts=ts) # + northward/ eastward surface shear stresses
    # Humidity replace at 1000hPa level results .
end

## Example (Assuming we have data for a given location)
## Populate state variables from ERA5 information. 

@info "Example of Unstable LES Boundary Layer"
z_int = FT(data_1.z)
z_sfc = FT(0)
z0m, z0b = FT(1.0f-2), FT(1.0f-3)
ts_int = data_1.ts
ts_sfc = data_sfc.ts
sc = SF.ValuesOnly{FT}(;
                       state_in = SF.InteriorValues(z_int, (FT(data_1.u), FT(data_1.v)), ts_int),
                       state_sfc = SF.SurfaceValues(z_sfc, (FT(0), FT(0)), ts_sfc),
                       z0m = z0m,
                       z0b = z0b,
                       )
@info "Evaluating Surface Fluxes using MOST"
result_fd = SF.surface_conditions(surface_params, sc, SF.FVScheme(); noniterative_stable_sol = false)
result_fv = SF.surface_conditions(surface_params, sc, SF.FDScheme(); noniterative_stable_sol = false)
#@info "UBL-Nishizawa, FD:", result_fd.L_MO, result_fd.shf, result.ρτxz
#@info "UBL-Nishizawa, FV:", result_fv.L_MO, result_fv.shf, result.ρτxz