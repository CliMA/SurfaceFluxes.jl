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
c_rsl = FT(0.7)         # RSL strength (Bonan supplemental programs use 0.7)

# φ(ζ) with Businger–Dyer (same universal functions as the solver)
phi(δz, transport) = UF.phi(uf_params, δz / L_MO, transport)

# RSL enhancement μ(δz): 1 above the RSL; model-specific form within.
mu_most(δz) = one(FT)
function mu_pg(δz)  # Physick & Garratt (1995) linear ramp
    ifelse(δz >= z_RSL, one(FT), one(FT) - c_rsl * (one(FT) - δz / z_RSL))
end
function mu_hf(δz)  # Harman & Finnigan (2007) / Bonan exponential
    ifelse(δz >= z_RSL, one(FT), exp(-c_rsl * (one(FT) - δz / z_RSL)))
end

# Integrate (φ μ / z) / κ from z0 to Δz with a composite midpoint rule.
function integrate_dimless(Δz, z0, transport, mu)
    z0_safe = max(z0, eps(FT))
    Δz_safe = max(Δz, z0_safe)
    n = 256
    # log-spaced nodes resolve the near-surface 1/z singularity cleanly
    ζ_nodes = range(log(z0_safe), log(Δz_safe); length = n + 1)
    acc = zero(FT)
    for i in 1:n
        z_lo = exp(ζ_nodes[i])
        z_hi = exp(ζ_nodes[i + 1])
        z_mid = sqrt(z_lo * z_hi)           # geometric midpoint
        acc += phi(z_mid, transport) * mu(z_mid) * log(z_hi / z_lo)
    end
    return acc / κ
end

# Profile from the roughness length. RSL curves are shifted so they coincide
# with MOST for z ≥ z_* (as drawn in Bonan Fig. 6.8).
function profile(z0, zs, transport, mu; match_at_rsl_top = false)
    vals = map(z -> integrate_dimless(z - d, z0, transport, mu), zs)
    if match_at_rsl_top
        most_star = integrate_dimless(z_star - d, z0, transport, mu_most)
        rsl_star = integrate_dimless(z_star - d, z0, transport, mu)
        vals = vals .+ (most_star - rsl_star)
    end
    return vals
end

zs_u = collect(range(d + z0m, FT(50); length = 200))
zs_θ = collect(range(d + z0h, FT(50); length = 200))

u_most = profile(z0m, zs_u, UF.MomentumTransport(), mu_most)
u_pg = profile(z0m, zs_u, UF.MomentumTransport(), mu_pg; match_at_rsl_top = true)
u_hf = profile(z0m, zs_u, UF.MomentumTransport(), mu_hf; match_at_rsl_top = true)

θ_most = profile(z0h, zs_θ, UF.HeatTransport(), mu_most)
θ_pg = profile(z0h, zs_θ, UF.HeatTransport(), mu_pg; match_at_rsl_top = true)
θ_hf = profile(z0h, zs_θ, UF.HeatTransport(), mu_hf; match_at_rsl_top = true)

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
    u_pg,
    zs_u;
    lw = 2,
    color = :dodgerblue,
    ls = :dash,
    label = "Physick–Garratt",
)
Plots.plot!(
    p1,
    u_hf,
    zs_u;
    lw = 2,
    color = :crimson,
    ls = :dot,
    label = "Harman–Finnigan",
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
    θ_pg,
    zs_θ;
    lw = 2,
    color = :dodgerblue,
    ls = :dash,
    label = "Physick–Garratt",
)
Plots.plot!(
    p2,
    θ_hf,
    zs_θ;
    lw = 2,
    color = :crimson,
    ls = :dot,
    label = "Harman–Finnigan",
)

rsl_fig = Plots.plot(p1, p2; layout = (1, 2), size = (900, 480))
Plots.savefig("RSL_profiles.svg")
nothing
