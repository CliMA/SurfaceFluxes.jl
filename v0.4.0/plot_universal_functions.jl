import SurfaceFluxes
const SF = SurfaceFluxes
const SFP = SF.Parameters
const UF = SF.UniversalFunctions

import CLIMAParameters
const CP = CLIMAParameters

include(joinpath(pkgdir(SurfaceFluxes), "parameters", "create_parameters.jl"))
const FT = Float32;
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")

import Plots


ζ = FT(-0.1):FT(0.001):FT(0.1);

L = FT(10);

ufts = (UF.GryanikType(), UF.BusingerType(), UF.GrachevType())

universal_functions(uft) = UF.universal_func(uft, L, create_uf_parameters(toml_dict, uft))

function save_ϕ_figs(ufts, ζ; ylims = nothing, fig_prefix = "", xaxis = :identity, yaxis = :identity)
    Plots.plot()
    for uft in ufts
        uf = universal_functions(uft)
        ϕ_m = UF.phi.(uf, ζ, UF.MomentumTransport())
        label = "$(typeof(uf).name)"
        Plots.plot!(ζ, ϕ_m; xlabel = "ζ", ylabel = "ϕ_m", label, ylims, xaxis, yaxis)
    end
    Plots.savefig("$(fig_prefix)_phi_m.svg")
    Plots.plot()
    for uft in ufts
        uf = universal_functions(uft)
        ϕ_h = UF.phi.(uf, ζ, UF.HeatTransport())
        label = "$(typeof(uf).name)"
        Plots.plot!(ζ, ϕ_h; xlabel = "ζ", ylabel = "ϕ_h", label, ylims, xaxis, yaxis)
    end
    Plots.savefig("$(fig_prefix)_phi_h.svg")
end
function save_ψ_figs(ufts, ζ; ylims = nothing, fig_prefix = "", xaxis = :identity, yaxis = :identity)
    Plots.plot()
    for uft in ufts
        uf = universal_functions(uft)
        ψ_m = UF.psi.(uf, ζ, UF.MomentumTransport())
        label = "$(typeof(uf).name)"
        Plots.plot!(ζ, ψ_m; xlabel = "ζ", ylabel = "ψ_m", label, ylims, xaxis, yaxis)
    end
    Plots.savefig("$(fig_prefix)_psi_m.svg")
    Plots.plot()
    for uft in ufts
        uf = universal_functions(uft)
        ψ_h = UF.psi.(uf, ζ, UF.HeatTransport())
        label = "$(typeof(uf).name)"
        Plots.plot!(ζ, ψ_h; xlabel = "ζ", ylabel = "ψ_h", label, ylims, xaxis, yaxis)
    end
    Plots.savefig("$(fig_prefix)_psi_h.svg")
end


save_ϕ_figs(ufts, FT(0):FT(0.01):FT(15); ylims = (0, 30), fig_prefix = "Gryanik12")
save_ψ_figs(ufts, FT(0):FT(0.01):FT(15); ylims = (-25, 0), fig_prefix = "Gryanik12")


save_ϕ_figs(
    ufts,
    10 .^ (FT(-3):0.1:FT(2));
    ylims = (0.1, 10^2),
    xaxis = :log10,
    yaxis = :log10,
    fig_prefix = "Gryanik3",
)


save_ϕ_figs(ufts, FT(-2.5):FT(0.01):FT(2); ylims = (-1, 8), fig_prefix = "Businger")
