using SurfaceFluxes, Documenter
using DocumenterCitations

# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "nul"

bib = CitationBibliography(joinpath(@__DIR__, "bibliography.bib"))


# Generate plots
cd(joinpath(@__DIR__, "src")) do
    include(joinpath(@__DIR__, "src", "plot_universal_functions.jl"))
    include(joinpath(@__DIR__, "src", "plot_bonan_profiles.jl"))
end

#! format: off
pages = Any[
    "Home" => "index.md",
    "Surface Fluxes Theory" => "SurfaceFluxes.md",
    "Universal Functions" => "UniversalFunctions.md",
    "Physical Scales" => "PhysicalScales.md",
    "Exchange Fluxes" => "ExchangeFluxes.md",
    "Prescribed Conditions" => "PrescribedConditions.md",
    "Test Suite" => "TestSuite.md",
    "API Reference" => "API.md",
]

mathengine = MathJax(Dict(
    :TeX => Dict(
        :equationNumbers => Dict(:autoNumber => "AMS"),
        :Macros => Dict(),
    ),
))

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    mathengine = mathengine,
    collapselevel = 1,
)
#! format: on

makedocs(
    bib,
    sitename = "SurfaceFluxes.jl",
    strict = false,
    format = format,
    checkdocs = :exports,
    clean = true,
    doctest = true,
    modules = [SurfaceFluxes],
    pages = pages,
)

deploydocs(
    repo = "github.com/CliMA/SurfaceFluxes.jl.git",
    target = "build",
    push_preview = true,
    devbranch = "main",
    forcepush = true,
)
