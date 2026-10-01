if !("." in LOAD_PATH) # for ease of local testing
    push!(LOAD_PATH, ".")
end

import ClimaParams as CP
import SurfaceFluxes
import SurfaceFluxes.UniversalFunctions
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters.SurfaceFluxesParameters
import Plots

FT = Float32
param_set = SurfaceFluxesParameters(FT, UniversalFunctions.BusingerParams)
uf_params = param_set.ufp
κ = param_set.von_karman_const

# --- Bonan (2019) Fig. 6.8 setup (same roughness as Fig. 6.4) -------------------
# Calculations as in Fig. 6.4, but with L_MO = -20 m and z_* = 49 m.
d = FT(19)              # displacement height [m]
z0m = FT(0.6)           # momentum roughness length [m]
z0h = FT(0.0816)        # scalar roughness length [m]
L_MO = FT(-20)          # Monin-Obukhov length [m]
z_star = FT(49)         # height of RSL top above ground [m]
z_RSL = z_star - d      # RSL depth above displacement height [m]
c_rsl = FT(0.7)         # RSL strength (Garratt 1980; Physick & Garratt 1995)

# RSL models anchored at the RSL top, so that the profiles coincide with MOST for
# z ≥ z_* (as drawn in Bonan Fig. 6.8). The exponential RSL factor with c = 0.7 is the
# form of Garratt (1980) and Physick & Garratt (1995).
rsl_exp = SurfaceFluxes.ExponentialRSL(FT; c_m = c_rsl, c_h = c_rsl, z_RSL = z_RSL)
rsl_lin = SurfaceFluxes.LinearRSL(FT; c_m = c_rsl, c_h = c_rsl, z_RSL = z_RSL)
no_rsl = SurfaceFluxes.NoRoughnessSubLayer()

# Dimensionless profile F̂(z)/κ = u(z)/u⋆ (or (θ(z) − θₛ)/θ⋆) at heights zs above ground
profile(z0, zs, transport, rsl) =
    map(zs) do z
        Δz_eff = z - d
        SurfaceFluxes.rsl_corrected_profile(
            uf_params,
            rsl,
            Δz_eff,
            Δz_eff / L_MO,
            z0,
            transport,
            UF.PointValueScheme(),
        ) / κ
    end

zs_u = collect(range(d + FT(1.05) * z0m, FT(50); length = 200))
zs_θ = collect(range(d + FT(1.05) * z0h, FT(50); length = 200))

u_most = profile(z0m, zs_u, UF.MomentumTransport(), no_rsl)
u_lin = profile(z0m, zs_u, UF.MomentumTransport(), rsl_lin)
u_exp = profile(z0m, zs_u, UF.MomentumTransport(), rsl_exp)

θ_most = profile(z0h, zs_θ, UF.HeatTransport(), no_rsl)
θ_lin = profile(z0h, zs_θ, UF.HeatTransport(), rsl_lin)
θ_exp = profile(z0h, zs_θ, UF.HeatTransport(), rsl_exp)

margin_kw = (; bottom_margin = 12Plots.mm, left_margin = 8Plots.mm, top_margin = 4Plots.mm)
axis_kw = (;
    framestyle = :box,
    tick_direction = :in,
    grid = false,
    legend = :topleft,
    legendfontsize = 8,
    ylim = (15, 50),
    margin_kw...,
)

# --- (a) Wind speed -------------------------------------------------------------
p1 = Plots.plot(
    u_most,
    zs_u;
    lw = 2,
    color = :black,
    ls = :solid,
    label = "MOST (no RSL)",
    xlim = (0, 7),
    xticks = 0:7,
    yticks = 15:5:50,
    xlabel = "u(z)/u⋆",
    ylabel = "Height (m)",
    title = "(a)",
    titlelocation = :left,
    axis_kw...,
)
Plots.plot!(
    p1,
    u_lin,
    zs_u;
    lw = 2,
    color = :dodgerblue,
    ls = :dash,
    label = "Linear RSL",
)
Plots.plot!(
    p1,
    u_exp,
    zs_u;
    lw = 2,
    color = :crimson,
    ls = :dot,
    label = "Exponential RSL (Garratt)",
)

# --- (b) Potential temperature --------------------------------------------------
p2 = Plots.plot(
    θ_most,
    zs_θ;
    lw = 2,
    color = :black,
    ls = :solid,
    label = "MOST (no RSL)",
    xlim = (0, 10),
    xticks = 0:2:10,
    yticks = 15:5:50,
    xlabel = "[θ(z) − θₛ]/θ⋆",
    ylabel = "",
    title = "(b)",
    titlelocation = :left,
    axis_kw...,
)
Plots.plot!(
    p2,
    θ_lin,
    zs_θ;
    lw = 2,
    color = :dodgerblue,
    ls = :dash,
    label = "Linear RSL",
)
Plots.plot!(
    p2,
    θ_exp,
    zs_θ;
    lw = 2,
    color = :crimson,
    ls = :dot,
    label = "Exponential RSL (Garratt)",
)

rsl_fig = Plots.plot(p1, p2; layout = (1, 2), size = (900, 480))
Plots.savefig("RSL_profiles.svg")
nothing
