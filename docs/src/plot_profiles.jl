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

""" 
Here, we reproduce Figure 6.4 from Bonan (2019) Chapter 6.

# References
- [Bonan2019](@cite) (Chapter 6 Figure 6.4)

# Original Research
- [Physick2019](@cite)

"""
T₀ = FT(273.15)
z_star = FT(49)
d = FT(19)
h_c = FT(22)
z0m = FT(0.6) # Figure 6.4
z0b = FT(0.135) * z0m

u_star_stable = FT(0.13)
u_star_unstable = FT(0.4)
θ_star_stable = FT(0.06)
θ_star_unstable = FT(-0.5)

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
Z = collect(range(FT(19.1), stop = FT(100), length = 250))

"""
save_profile(param_set, sc, ca, L_MOs, Z, X_in, X_sfc, transport, uft, scheme, x_star, d)

Saves profiles of variable X given values of Z coordinates. Follows Nishizawa equation (21,22)

## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - L_MOs: Monin-Obukhov length(s)
  - ca: Canopy Type (e.g. NoCanopy, SparseCanopy, DenseCanopy)
  - Z: Z coordinate(s) (within surface layer) for which variable values are required
  - X_in, X_sfc: For variable X, values at interior and surface nodes
  - transport: Transport type, (e.g. Momentum or Heat, used to determine physical scale coefficients)
  - uft: A Universal Function type, (returned by, e.g., Businger())
  - scheme: Discretization scheme (currently supports FD and FV)
  - x_star: characteristic scale for variable x
  - d: Displacement height (measure of the spatial lengthscale of the effect of the canopy roughness on near-wall turbulence)

"""
function save_profile(
    param_set::SF.APS,
    sc::SF.AbstractSurfaceConditions,
    ca::Union{SF.NoCanopy, SF.SparseCanopy, SF.DenseCanopy},
    L_MOs::Array{FT, 1},
    Z::Array{FT, 1},
    X_in::FT,
    X_sfc::FT,
    transport,
    uft::UF.AUFT,
    scheme::Union{SF.FVScheme, SF.FDScheme},
    x_star::FT,
    d::FT;
    xlims = nothing,
    ylims = nothing,
    xlabel = L"$\frac{u}{u_{\star}}$",
    ylabel = L"$z$",
    fig_prefix = "",
    xaxis = :identity,
    yaxis = :identity,
)
    Plots.plot()
    for L_MO in L_MOs
        x_i = map(Z) do z
            Zi = typeof(ca) == SF.NoCanopy ? FT(z - d) : FT(z)
            dx = SF.recover_profile(param_set, sc, ca, L_MO, Zi, X_in, X_sfc, transport, uft, scheme)
        end

        Δx = κ * (x_i .- X_sfc) ./ x_star

        Plots.plot!(Δx, Z, label = L"L_{MO} = %$L_MO")
        Plots.plot!(; xlabel, ylabel, ylims, xlims, grid = :off, legend = :outerright)

        Plots.savefig("$(fig_prefix)_profile.svg")
    end
end



# Save profile for unstable and stable L_MO set
Stable_L_MOs = FT[30, 50, 1000]
Unstable_L_MOs = FT[-10, -50, -1000]
testcanopy = SurfaceFluxes.SparseCanopy{FT}(d, z_star)
nocanopy = SurfaceFluxes.NoCanopy()

# Bonan2019 Fig. 6.4a (With Canopy)
save_profile(
    param_set,
    sc,
    testcanopy,
    Unstable_L_MOs,
    Z,
    FT(3.7),
    FT(0),
    UF.MomentumTransport(),
    uft,
    SF.FDScheme(),
    u_star_unstable,
    d;
    xlims = (0, 4),
    ylims = (15, 50),
    fig_prefix = "Fig6.4a_canopy",
)

# Bonan2019 Fig. 6.4a (No Canopy)
save_profile(
    param_set,
    sc,
    nocanopy,
    Unstable_L_MOs,
    Z,
    FT(3.7),
    FT(0),
    UF.MomentumTransport(),
    uft,
    SF.FDScheme(),
    u_star_unstable,
    d;
    xlims = (0, 4),
    ylims = (15, 50),
    fig_prefix = "Fig6.4a",
)

# Bonan2019 Fig. 6.4b (With Canopy)
save_profile(
    param_set,
    sc,
    testcanopy,
    Unstable_L_MOs,
    Z,
    FT(298.0),
    θ_sfc,
    UF.HeatTransport(),
    uft,
    SF.FDScheme(),
    θ_star_unstable,
    d;
    xlims = (-8, 0),
    ylims = (15, 50),
    xlabel = L"$\frac{\theta}{\theta_{\star}}$",
    fig_prefix = "Fig6.4b_canopy",
)

# Bonan2019 Fig. 6.4b (No Canopy)
save_profile(
    param_set,
    sc,
    nocanopy,
    Unstable_L_MOs,
    Z,
    FT(298.0),
    θ_sfc,
    UF.HeatTransport(),
    uft,
    SF.FDScheme(),
    θ_star_unstable,
    d;
    xlims = (-8, 0),
    ylims = (15, 50),
    xlabel = L"$\frac{\theta}{\theta_{\star}}$",
    fig_prefix = "Fig6.4b",
)

# Bonan2019 Fig. 6.4c (With Canopy)
save_profile(
    param_set,
    sc,
    testcanopy,
    Stable_L_MOs,
    Z,
    FT(3.7),
    FT(0),
    UF.MomentumTransport(),
    uft,
    SF.FDScheme(),
    u_star_stable,
    d;
    xlims = (0, 4),
    ylims = (15, 50),
    fig_prefix = "Fig6.4c_canopy",
)

# Bonan2019 Fig. 6.4c (No Canopy)
save_profile(
    param_set,
    sc,
    nocanopy,
    Stable_L_MOs,
    Z,
    FT(3.7),
    FT(0),
    UF.MomentumTransport(),
    uft,
    SF.FDScheme(),
    u_star_stable,
    d;
    xlims = (0, 4),
    ylims = (15, 50),
    fig_prefix = "Fig6.4c",
)

# Bonan2019 Fig. 6.4d (With Canopy)
save_profile(
    param_set,
    sc,
    testcanopy,
    Stable_L_MOs,
    Z,
    FT(298.0),
    θ_sfc,
    UF.HeatTransport(),
    uft,
    SF.FDScheme(),
    θ_star_stable,
    d;
    xlims = (0, 2),
    ylims = (15, 50),
    xlabel = L"$\frac{\theta}{\theta_{\star}}$",
    fig_prefix = "Fig6.4d_canopy",
)

# Bonan2019 Fig. 6.4d (No Canopy)
save_profile(
    param_set,
    sc,
    nocanopy,
    Stable_L_MOs,
    Z,
    FT(298.0),
    θ_sfc,
    UF.HeatTransport(),
    uft,
    SF.FDScheme(),
    θ_star_stable,
    d;
    xlims = (0, 2),
    ylims = (15, 50),
    xlabel = L"$\frac{\theta}{\theta_{\star}}$",
    fig_prefix = "Fig6.4d",
)


# tol_neutral = FT(SFP.cp_d(param_set) / 100);
# z_level = FT(5);
# ζ = collect(range(-5, stop = 1, length = 50));
# CD = [];
# for i in 1:length(ζ)
#     push!(CD, SF.momentum_exchange_coefficient(param_set, ζ[i] * z_level, sc, uft, SF.FDScheme(), tol_neutral))
# end
