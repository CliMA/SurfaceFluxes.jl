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
FT = Float64
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
surface_params = create_parameters(toml_dict, UF.Businger())
thermo_params = SFP.thermodynamics_params(surface_params)

# 950hPa level
extract_var(ds, n; 
            ilon= 1, 
            ilat= 1, 
            it= 1) = ds[n][:,:,it] # * ds[n].attrib["scale_factor"] + ds[n].attrib["add_offset"] # NCDatasets deals with this

# Check in Loop
data_int = NCDataset("era5_97500.nc","r") do ds
  T = extract_var(ds, "t") # temperature
  u = extract_var(ds, "u") # eastward velo
  v = extract_var(ds, "v") # northward velo
  z = extract_var(ds, "z")/ 9.81 # altitude from "z" geopotential
  q = extract_var(ds, "q") # total moisture (specific humidity)
  p = 97500.0 # source data pressure level (950hPa)
  ts = @. TD.PhaseEquil_pTq(thermo_params, p, T, q)
  ρ = @. air_density(thermo_params, ts)
  (; T = T, u = u, v = v, z = z, q = q, ρ = ρ, 
  ts=ts)
end

# 2-m / surface level
data_sfc = NCDataset("era5_surface.nc","r") do ds
    T = extract_var(ds, "t2m") # temperature
    p = extract_var(ds, "sp")
    u = 0.0 # eastward velo (may be non zero for ocean)
    v = 0.0 # northward velo
    z = 0.0 # altitude from "z" geopotential
    e = @. saturation_vapor_pressure(thermo_params, T, TD.Liquid())
    # tetens_eq(𝒯) = 0.61078 * exp(17.27 * 𝒯  / (𝒯 + 237.3)) 
    # e = tetens_eq(T)
    # TODO: Update with ClimaParameters values 
    # (EarthParameters <: AbstractParameterSet)
    R = 287.1
    R_m = 461.5
    #q_sat = R / R_m * e * (p - ((1 - R / R_m) * e)) # saturated total moisture (specific humidity) # BUG
    
    # q_sat following ECMWF documentation (TODO Replace with CLIMAParameters)
    Rdry = 287.0597 
    Rvap=461.5250 
    a₁=611.21
    a₃=17.502
    a₄=32.19
    T₀=273.16
    sat_vap  = @. a₁ * exp(a₃ * (T - T₀) / (T - a₄))
    q_sat= @. (Rdry/Rvap) * sat_vap / (p - ((1 - Rdry / Rvap) * sat_vap))

    ts = @. TD.PhaseEquil_pTq(thermo_params, p,T, q_sat)
    ρ = @. air_density(thermo_params, ts)

    (;shf=extract_var(ds, "sshf"), 
      lhf=extract_var(ds, "slhf"), 
      evap = extract_var(ds, "e"), 
      ustar = extract_var(ds,"zust"),
      nsss = extract_var(ds, "nsss"),
      ewss = extract_var(ds, "ewss"),
      T = T, u = u, v = v, z = z, q = q_sat, ρ = ρ, 
      ts=ts) 
    # Humidity replace at 1000hPa level results .
end

for (ii, ilon) in enumerate((1))
  for (jj, ilat) in enumerate((1))
    ## Example (Assuming we have data for a given location)
    ## Populate state variables from ERA5 information. 
    z_int = FT(data_int.z[ilon, ilat])
    z_sfc = FT(0)
      z0m, z0b = FT(1.0f-2), FT(1.0f-3)
      ts_int = data_int.ts[ilon, ilat]
      ts_sfc = data_sfc.ts[ilon, ilat]
      sc = SF.ValuesOnly{FT}(;
              state_in = SF.InteriorValues(z_int, (FT(data_int.u[ilon, ilat]), FT(data_int.v[ilon,ilat])), ts_int),
              state_sfc = SF.SurfaceValues(z_sfc, (FT(0), FT(0)), ts_sfc),
              z0m = z0m,
              z0b = z0b,
           )
       result_fd = SF.surface_conditions(surface_params, sc, SF.FVScheme(); noniterative_stable_sol = false);
  #    result_fv = SF.surface_conditions(surface_params, sc, SF.FDScheme(); noniterative_stable_sol = false);
  #    @info "FD:", result_fd.shf, result_fd.lhf, result_fd.ρτxz, result_fd.ρτyz, result_fd.ustar
      @info "FD:"  result_fd.ustar
  #    @info "RAW:", data_sfc.shf/3600, data_sfc.lhf/3600, data_sfc.ewss/3600, data_sfc.nsss/3600, data_sfc.ustar
      @info "RAW:" data_sfc.ustar[ilon, ilat]
  end
end
      
