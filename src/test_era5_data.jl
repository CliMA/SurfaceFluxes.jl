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
using Plots

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
# TODO: Update all constants with corresponding climaparameters values

extract_var(ds, n; 
            ilon= 1, 
            ilat= 1, 
            it= 1) = ds[n][:,:,it] # * ds[n].attrib["scale_factor"] + ds[n].attrib["add_offset"] 
# NCDatasets determines the scaling + offsets by default

# Check in Loop
data_int = NCDataset("era5_97500.nc","r") do ds
  T = extract_var(ds, "t") # temperature
  u = extract_var(ds, "u") # eastward velo
  v = extract_var(ds, "v") # northward velo
  z = extract_var(ds, "z")/ 9.81 # altitude from ERA5 "z" geopotential
  q = extract_var(ds, "q") # total moisture (specific humidity)
  p = 97500.0 # source data pressure level (975hPa)
  ts = @. TD.PhaseEquil_pTq(thermo_params, p, T, q) # determine thermodynamic state (assumes no liquid-ice content)
  ρ = @. air_density(thermo_params, ts) # compute density from thermo-state
  (; T = T, u = u, v = v, z = z, q = q, ρ = ρ, ts=ts) # return variables
end

# 2-m / surface level
data_sfc = NCDataset("era5_surface.nc","r") do ds
    T = extract_var(ds, "t2m") # temperature
    p = extract_var(ds, "sp")
    u = 0.0 # eastward velo (may be non zero for ocean)
    v = 0.0 # northward velo
    e = @. saturation_vapor_pressure(thermo_params, T, TD.Liquid())
    R = 287.1
    R_m = 461.5
    # (Replace with CLIMAParameters)
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
      z = extract_var(ds, "z")/ 9.81, # altitude from ERA5 "z" geopotential
      z0m = extract_var(ds, "fsr"), # Raw data stores the estimated aerodynamic roughness length for momentum!
      z0b = exp.(extract_var(ds, "flsr")), # Raw data stores natural log of the roughness length for heat!
      λ = ds["longitude"][:], 
      𝜃 = ds["latitude"][:], 
      T = T, u = u, v = v, q = q_sat, ρ = ρ, 
      ts=ts) 
end

ustar = zeros(1440,721)
shf = zeros(1440,721)
lhf = zeros(1440,721)
evaporation = zeros(1440,721)
altitude = zeros(1440,721)
lmo = zeros(1440,721)
ρτxz = zeros(1440,721)
ρτyz = zeros(1440,721)
for (ii, ilon) in enumerate(data_sfc.λ)
  for (jj, ilat) in enumerate(data_sfc.𝜃)
    ## Example (Assuming we have data for a given location)
    ## Populate state variables from ERA5 information. 
    z_sfc = FT(data_sfc.z[ii, jj])
    Δz_grid = FT(100)
    z_int = z_sfc .+ Δz_grid # FT(data_int.z[ii, jj])
    # ERA fields (local values of fsr and flsr)
    z0m = FT(0.001)
    z0b = FT(0.001)
    ts_int = data_int.ts[ii, jj]
    ts_sfc = data_sfc.ts[ii, jj]
    sc = SF.ValuesOnly{FT}(;
            state_in = SF.InteriorValues(z_int, (FT(data_int.u[ii, jj]), FT(data_int.v[ii,jj])), ts_int),
            state_sfc = SF.SurfaceValues(z_sfc, (FT(0), FT(0)), ts_sfc),
            z0m = z0m,
            z0b = z0b,
         )
    result_fd = SF.surface_conditions(surface_params, sc, SF.FVScheme(); noniterative_stable_sol = false);
    result_fv = SF.surface_conditions(surface_params, sc, SF.FDScheme(); noniterative_stable_sol = false);
    ustar[ii,jj] = result_fd.ustar 
    lmo[ii,jj] = result_fd.L_MO # Monin-Obukhov
    shf[ii,jj] = result_fd.shf
    lhf[ii,jj] = result_fd.lhf
    evaporation[ii,jj] = result_fd.evaporation
    ρτxz[ii,jj] = result_fd.ρτxz # Eastward turbulent shear stresses
    ρτyz[ii,jj] = result_fd.ρτyz # Northward turbulene shear stresses
  end
end

Plots.heatmap(reverse(ustar[:,:,1]', dims=1), clims=(-0, 1.5), c = :RdGy_9)
