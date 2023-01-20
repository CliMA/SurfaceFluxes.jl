using Plots
using LaTeXStrings

import SurfaceFluxes
SurfaceFluxes.error_on_non_convergence() = true
const SF = SurfaceFluxes
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

### Test Unstablea
# The material in this section follows Chapter 6 in Bonan (2019). 
# The source material is Physick and Garratt (1995).
T₀ = FT(273.15)
L_MO = FT(-20)
z_star = FT(49)
d = FT(19)
h_c  = FT(22)
z0m = FT(0.6) # Figure 6.4
z0b = FT(0.135) * z0m
u_star = FT(0.6)
θ_star = FT(-0.5)
κ = 0.4

# For these measured values, see Chapter 6.6 page 89 of Bonan (2019)
z_measured = FT[21, 29, 21];
u_measured = FT[1.0, 2.1, 1.0];
θ_measured = FT[29.0, 28.1, 29.5] .+ T₀

p_sfc = FT(99340) # Pa
q_sfc = FT(0.0107) # kg/kg

θ_sfc = FT(289.7) # K
ts_sfc_test = Thermodynamics.PhaseEquil_pθq(thermo_params, 99340.0, 289.7, 0.0107)
ts_int_test = Thermodynamics.PhaseEquil_pθq(thermo_params, 95342.0, 298.0, 0.0085)
state_in = SF.InteriorValues(FT(350), (FT(0.0), FT(0)), ts_int_test)
state_sfc = SF.SurfaceValues(FT(0), (FT(0), FT(0)), ts_sfc_test)
sc = SF.ValuesOnly{FT}(; state_in, state_sfc, z0m, z0b)
rsc = SF.surface_conditions(param_set, sc, SF.FDScheme())
Z = collect(range(FT(25), stop=FT(100), length=100))

u_i_canopy = [];
θ_i_canopy = [];

testcanopy = SurfaceFluxes.SparseCanopy{FT}(d,z_star)
nocanopy = SurfaceFluxes.NoCanopy()
for (iz, z) in enumerate(Z)
  du = SF.recover_profile(param_set, 
                          sc, 
                          testcanopy,
                          L_MO, 
                          FT(z), 
                          FT(3.7),
                          FT(0),
                          UF.MomentumTransport(), 
                          uft, 
                          SF.FDScheme())
  dθ = SF.recover_profile(param_set, 
                            sc, 
                            testcanopy,
                            L_MO, 
                            FT(z), 
                            FT(298.0), 
                            θ_sfc,
                            UF.HeatTransport(), 
                            uft, 
                            SF.FDScheme())
  push!(u_i_canopy, du)
  push!(θ_i_canopy, dθ)
end

u_i= [];
θ_i= [];
for (iz, z) in enumerate(Z)
  du = SF.recover_profile(param_set, 
                          sc, 
                          nocanopy,
                          L_MO, 
                          FT(z - d), 
                          FT(3.7),
                          FT(0),
                          UF.MomentumTransport(), 
                          uft, 
                          SF.FDScheme())
  push!(u_i, du)
  dθ = SF.recover_profile(param_set, 
                            sc, 
                            nocanopy,
                            L_MO, 
                            FT(z - d), 
                            FT(298.0),
                            θ_sfc,
                            UF.HeatTransport(), 
                            uft, 
                            SF.FDScheme())
  push!(θ_i, dθ)
end
kΔu = κ * u_i ./ u_star;
kΔθ = κ * (θ_i .- θ_sfc) ./ θ_star;
kΔθ_canopy = κ * (θ_i_canopy .- θ_sfc) ./ θ_star;

Plots.plot(u_i / u_star,Z, m = :sq, label="No-RSL")
Plots.plot!(u_i_canopy ./ u_star,Z, m = :o, label="RSL")
Plots.plot!(;
            xlabel = L"$\frac{u}{u_{\star}}$",
            ylabel = L"$z$",
            ylim=(0, 60), 
            grid=:off,
            legend=:outerright)

#Plots.plot(kΔθ / θ_star,Z, m = :sq, label="No-RSL")
#Plots.plot!(kΔθ_canopy ,Z, m = :o, label="RSL")
#Plots.plot!(;
#            xlabel = L"$\frac{\theta}{\theta_{\star}}$",
#            ylabel = L"$z$",
#            ylim=(0, 60), 
#            grid=:off,
#            legend=:outerright)


tol_neutral = FT(SFP.cp_d(param_set) / 100);
z_level = FT(5);
ζ = collect(range(-5,stop=1, length=50));
CD = [];
for i in 1:length(ζ)
  push!(CD, SF.momentum_exchange_coefficient(
    param_set, 
    ζ[i] * z_level,
    sc,
    uft,
    SF.FDScheme(),
    tol_neutral,
   ))
end
