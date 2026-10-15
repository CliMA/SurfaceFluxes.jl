using Plots

using SurfaceFluxes
const SF = SurfaceFluxes
SurfaceFluxes.error_on_non_convergence() = true

import SurfaceFluxes.UniversalFunctions as UF
import Thermodynamics
Thermodynamics.print_warning() = false

include(joinpath(pkgdir(SurfaceFluxes), "parameters", "create_parameters.jl"))
const SFP = SurfaceFluxes.Parameters
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

κ = SFP.von_karman_const(param_set)

# For these measured values, see Chapter 6.6 page 89 of Bonan (2019)
z_measured = FT[21, 29, 21];
u_measured = FT[1.0, 2.1, 1.0];
θ_measured = FT[29.0, 28.1, 29.5] .+ T₀

p_sfc = FT(99340) # Pa
q_sfc = FT(0.0107) # kg/kg

θ_sfc = FT(289.7) # K
ts_sfc_test = Thermodynamics.PhaseEquil_pθq(thermo_params, 99340.0, 289.7, 0.0107)
ts_int_test = Thermodynamics.PhaseEquil_pθq(thermo_params, 95342.0, 298.0, 0.0085)
state_in = SurfaceFluxes.InteriorValues(FT(100), (FT(1.0), FT(0)), ts_int_test)
state_sfc = SurfaceFluxes.SurfaceValues(FT(0), (FT(0), FT(0)), ts_sfc_test)
sc = SurfaceFluxes.ValuesOnly{FT}(; state_in, state_sfc, z0m, z0b)
ts_sfc_test = Thermodynamics.PhaseEquil_pθq(thermo_params, 100000.0, 289.7, 0.0)
ts_int_test = Thermodynamics.PhaseEquil_pθq(thermo_params, 99990.0, 298.0, 0.0)
Z = collect(range(FT(d + 10^-10), stop = FT(100), length = 250))

"""
save_profile(param_set, sc, ca, L_MOs, Z, X_sfc, transport, uft, scheme, x_star, d)

Saves profiles of variable X given values of Z coordinates. Follows Nishizawa equation (21,22)

## Arguments
  - param_set: Abstract Parameter Set containing physical, thermodynamic parameters.
  - sc: Container for surface conditions based on known combination
        of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment
  - L_MOs: Monin-Obukhov length(s)
  - Z: Z coordinate(s) (within surface layer) for which variable values are required
  - X_sfc: For variable X, values at interior and surface nodes
  - transport: Transport type, (e.g. Momentum or Heat, used to determine physical scale coefficients)
  - uft: A Universal Function type, (returned by, e.g., Businger())
  - scheme: Discretization scheme (currently supports FD and FV)
  - rsl : Roughness Sublayer Formulation (e.g. NoRSL, PhysickRSL, DeRidderRSL)
  - x_star: characteristic scale for variable x
  - d: Displacement height (measure of the spatial lengthscale of the effect of the canopy roughness on near-wall turbulence)
"""
function save_profile(
    param_set::SurfaceFluxes.APS,
    sc::SurfaceFluxes.AbstractSurfaceConditions,
    L_MOs::Array{FT, 1},
    Z::Array{FT, 1},
    X_sfc,
    transport,
    uft::UF.AUFT,
    scheme::Union{SurfaceFluxes.FVScheme, SurfaceFluxes.FDScheme},
    rsl::SurfaceFluxes.AbstractRoughnessSublayerType,
    x_star,
    d;
    title = nothing,
    xlims = nothing,
    ylims = nothing,
    xlabel = "u(z)",
    ylabel = "z",
    fig_prefix = "",
    xaxis = :identity,
    yaxis = :identity,
)
    Plots.plot()
    for L_MO in L_MOs
        x_i = map(Z) do z
            Zi = typeof(rsl) == SurfaceFluxes.NoRSL ? FT(z - d) : FT(z)
            state_in = SurfaceFluxes.InteriorValues(FT(Zi), (FT(1.0), FT(0)), ts_int_test)
            state_sfc = SurfaceFluxes.SurfaceValues(FT(0), (FT(0), FT(0)), ts_sfc_test)
            sc = SurfaceFluxes.ValuesOnly{FT}(; state_in, state_sfc, z0m, z0b)
            rsc = SurfaceFluxes.surface_conditions(param_set, sc, SurfaceFluxes.FDScheme())
            dx = SurfaceFluxes.recover_profile(param_set, sc, L_MO, Zi, X_sfc, x_star, transport, uft, scheme, rsl)
        end

        uf = UF.universal_func(uft, L_MO, SFP.uf_params(param_set))
        _π_group = FT(UF.π_group(uf, transport))

        Δx = @. (x_i - X_sfc)

        Plots.plot!(Δx, Z, label = "L_MO = $L_MO")
        Plots.plot!(; title, xlabel, ylabel, ylims, xlims, grid = :off, legend = :outerright, titlefontalign = :center)

        Plots.savefig("$(fig_prefix)_profile.png")
    end
end



# Save profile for unstable and stable L_MO set
Stable_L_MOs = FT[30, 50, 1000]
Unstable_L_MOs = FT[-10, -50, -1000]

testcanopy = SurfaceFluxes.SparseCanopy{FT}(d, z_star)
PhysickRSL = SurfaceFluxes.PhysickRSL(testcanopy)
NoRSL = SurfaceFluxes.NoRSL()


# Bonan2019 Fig. 6.4a (With PG95 RSL)
save_profile(
    param_set,
    sc,
    Unstable_L_MOs,
    Z,
    FT(0),
    UF.MomentumTransport(),
    uft,
    SurfaceFluxes.FDScheme(),
    PhysickRSL,
    u_star_unstable,
    d;
    xlims = (0, 4),
    ylims = (15, 50),
    fig_prefix = "Fig6.4a_canopy",
    title = "Vertical Profile of Wind Velocity (PG95 RSL)",
)

# Bonan2019 Fig. 6.4a (No RSL Model)
save_profile(
    param_set,
    sc,
    Unstable_L_MOs,
    Z,
    FT(0),
    UF.MomentumTransport(),
    uft,
    SurfaceFluxes.FDScheme(),
    NoRSL,
    u_star_unstable,
    d;
    xlims = (0, 4),
    ylims = (15, 50),
    fig_prefix = "Fig6.4a",
    title = "Vertical Profile of Wind Velocity (No RSL Model)",
)

# Bonan2019 Fig. 6.4b (With PG95 RSL)
save_profile(
    param_set,
    sc,
    Unstable_L_MOs,
    Z,
    θ_sfc,
    UF.HeatTransport(),
    uft,
    SurfaceFluxes.FDScheme(),
    PhysickRSL,
    θ_star_unstable,
    d;
    xlims = (-8, 0),
    ylims = (15, 50),
    xlabel = "θ - θ_sfc",
    fig_prefix = "Fig6.4b_canopy",
    title = "Vertical Profile of Temperature (PG95 RSL)",
)

# Bonan2019 Fig. 6.4b (No RSL Model)
save_profile(
    param_set,
    sc,
    Unstable_L_MOs,
    Z,
    θ_sfc,
    UF.HeatTransport(),
    uft,
    SurfaceFluxes.FDScheme(),
    NoRSL,
    θ_star_unstable,
    d;
    xlims = (-8, 0),
    ylims = (15, 50),
    xlabel = "θ- θ_sfc",
    fig_prefix = "Fig6.4b",
    title = "Vertical Profile of Temperature (No RSL Model)",
)

# Bonan2019 Fig. 6.4c (With PG95 RSL)
save_profile(
    param_set,
    sc,
    Stable_L_MOs,
    Z,
    FT(0),
    UF.MomentumTransport(),
    uft,
    SurfaceFluxes.FDScheme(),
    PhysickRSL,
    u_star_stable,
    d;
    xlims = (0, 4),
    ylims = (15, 50),
    fig_prefix = "Fig6.4c_canopy",
    title = "Vertical Profile of Wind Velocity (PG95 RSL)",
)

# Bonan2019 Fig. 6.4c (No RSL Model)
save_profile(
    param_set,
    sc,
    Stable_L_MOs,
    Z,
    FT(0),
    UF.MomentumTransport(),
    uft,
    SurfaceFluxes.FDScheme(),
    NoRSL,
    u_star_stable,
    d;
    xlims = (0, 4),
    ylims = (15, 50),
    fig_prefix = "Fig6.4c",
    title = "Vertical Profile of Wind Velocity (No RSL Model)",
)

# Bonan2019 Fig. 6.4d (With PG95 RSL)
save_profile(
    param_set,
    sc,
    Stable_L_MOs,
    Z,
    θ_sfc,
    UF.HeatTransport(),
    uft,
    SurfaceFluxes.FDScheme(),
    PhysickRSL,
    θ_star_stable,
    d;
    xlims = (0, 2),
    ylims = (15, 50),
    xlabel = "θ-θ_sfc",
    fig_prefix = "Fig6.4d_canopy",
    title = "Vertical Profile of Temperature (PG95 RSL)",
)

# Bonan2019 Fig. 6.4d (No RSL Model)
save_profile(
    param_set,
    sc,
    Stable_L_MOs,
    Z,
    θ_sfc,
    UF.HeatTransport(),
    uft,
    SurfaceFluxes.FDScheme(),
    NoRSL,
    θ_star_stable,
    d;
    xlims = (0, 2),
    ylims = (15, 50),
    xlabel = "θ-θ_sfc",
    fig_prefix = "Fig6.4d",
    title = "Vertical Profile of Temperature (No RSL Model)",
)
