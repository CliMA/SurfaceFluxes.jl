if !("." in LOAD_PATH) # for ease of local testing
    push!(LOAD_PATH, ".")
end

import Plots

import SurfaceFluxes as SF
import SurfaceFluxes.Parameters.SurfaceFluxesParameters as SFP
import SurfaceFluxes.UniversalFunctions as UF
import ClimaParams as CP

const FT = Float32

ζ = FT(-0.1):FT(0.001):FT(0.1)

L = FT(10)

ufps = (
    UF.GryanikParams(FT),
    UF.BusingerParams(FT),
    UF.GrachevParams(FT),
)

function save_ϕ_figs(
    ufps,
    ζ;
    ylims = nothing,
    fig_prefix = "",
    xaxis = :identity,
    yaxis = :identity,
)
    Plots.plot()
    for ufp in ufps
        ϕ_m = UF.phi.(ufp, ζ, UF.MomentumTransport())
        label = "$(typeof(ufp))"
        Plots.plot!(
            ζ,
            ϕ_m;
            xlabel = "ζ",
            ylabel = "ϕ_m",
            label,
            ylims,
            xaxis,
            yaxis,
        )
    end
    Plots.savefig("$(fig_prefix)_phi_m.svg")
    Plots.plot()
    for ufp in ufps
        ϕ_h = UF.phi.(ufp, ζ, UF.HeatTransport())
        label = "$(typeof(ufp))"
        Plots.plot!(
            ζ,
            ϕ_h;
            xlabel = "ζ",
            ylabel = "ϕ_h",
            label,
            ylims,
            xaxis,
            yaxis,
        )
    end
    Plots.savefig("$(fig_prefix)_phi_h.svg")
    return nothing
end
function save_ψ_figs(
    ufps,
    ζ;
    ylims = nothing,
    fig_prefix = "",
    xaxis = :identity,
    yaxis = :identity,
)
    Plots.plot()
    for ufp in ufps
        ψ_m = UF.psi.(ufp, ζ, UF.MomentumTransport())
        label = "$(typeof(ufp))"
        Plots.plot!(
            ζ,
            ψ_m;
            xlabel = "ζ",
            ylabel = "ψ_m",
            label,
            ylims,
            xaxis,
            yaxis,
        )
    end
    Plots.savefig("$(fig_prefix)_psi_m.svg")
    Plots.plot()
    for ufp in ufps
        ψ_h = UF.psi.(ufp, ζ, UF.HeatTransport())
        label = "$(typeof(ufp))"
        Plots.plot!(
            ζ,
            ψ_h;
            xlabel = "ζ",
            ylabel = "ψ_h",
            label,
            ylims,
            xaxis,
            yaxis,
        )
    end
    Plots.savefig("$(fig_prefix)_psi_h.svg")
    return nothing
end


# Gryanik Plots
save_ϕ_figs(
    ufps,
    FT(0):FT(0.01):FT(15);
    ylims = (0, 30),
    fig_prefix = "Gryanik12",
)
save_ψ_figs(
    ufps,
    FT(0):FT(0.01):FT(15);
    ylims = (-25, 0),
    fig_prefix = "Gryanik12",
)

save_ϕ_figs(
    ufps,
    10 .^ (FT(-3):0.1:FT(2));
    ylims = (0.1, 10^2),
    xaxis = :log10,
    yaxis = :log10,
    fig_prefix = "Gryanik3",
)

# Businger Plots
save_ϕ_figs(
    ufps,
    FT(-2.5):FT(0.01):FT(2);
    ylims = (-1, 8),
    fig_prefix = "Businger",
)

# Bonan Plots
save_ϕ_figs(ufps, FT(-2):FT(0.01):FT(1); ylims = (0, 6), fig_prefix = "Bonan")
save_ψ_figs(ufps, FT(-2):FT(0.01):FT(1); ylims = (-5, 4), fig_prefix = "Bonan")
