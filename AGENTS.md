# SurfaceFluxes.jl Agent Guide

## Ecosystem Guidelines

Please refer to the shared CliMA agent index for ecosystem-wide rules regarding architecture, performance, code quality, infrastructure, and workflows:

- [docs/dev-guides/AGENTS.md](docs/dev-guides/AGENTS.md) — Shared CliMA agent guidelines.

> Shared guides live at `docs/dev-guides/` and are vendored from the canonical source:
> <https://github.com/CliMA/DeveloperGuides>. Edit shared guides there, not here. They are
> synced automatically each month by `.github/workflows/update_dev_guides.yml`.

## Before You Act: Agent Autonomy

Before making changes that are externally visible or consequential (`git push`, version bumps, CI config changes, public API renames), check [docs/dev-guides/workflow/agent_autonomy.md](docs/dev-guides/workflow/agent_autonomy.md). The boundaries listed there require explicit user approval.

## Repo-Specific Guidelines

SurfaceFluxes.jl computes turbulent surface fluxes of momentum, heat, and moisture using Monin-Obukhov Similarity Theory (MOST). It is GPU-compatible, AD-compatible (ForwardDiff), and designed to be broadcast over arrays of heterogeneous surfaces.

### Architecture

- **Two modules**: the top-level `SurfaceFluxes` module (physics, solver, fluxes) and the nested `SurfaceFluxes.UniversalFunctions` submodule (stability functions ϕ, ψ, Ψ). `SurfaceFluxes.Parameters` holds the parameter struct and accessors.
- **Pure, stateless functions**: physics functions take a `param_set`, state, and configuration and return values or a [`SurfaceFluxConditions`](docs/src/API.md) struct; they do not mutate global state.
- **GPU/AD compatibility**: hot paths use `ifelse` rather than branches to avoid warp divergence; both branches of an `ifelse` must be evaluable (guard invalid inputs with `min`/`max` so the discarded branch never errors). Avoid keyword-argument constructors in kernels. Do not contaminate the `Float32` path with `Float64` literals — wrap constants in `FT(...)`. See [docs/dev-guides/performance/](docs/dev-guides/performance/).
- **Solver design**: `surface_fluxes` dispatches to one of four modes (prescribed coefficients, prescribed fluxes, prescribed heat+drag, or the iterative MOST solve). The iterative solve finds the stability parameter `ζ` as the root of `Ri_b(state) − Ri_b(ζ)` via [RootSolvers.jl](https://github.com/CliMA/RootSolvers.jl), using a fixed iteration count by default for branch-free GPU execution.

### Source layout

| Path | Purpose |
|------|---------|
| `src/SurfaceFluxes.jl` | Top-level module: `surface_fluxes` entry point, mode dispatch, the MOST solve (`solve_monin_obukhov`, `ResidualFunction`) |
| `src/UniversalFunctions.jl` | `UniversalFunctions` submodule: ϕ/ψ/Ψ for Businger, Gryanik, Grachev; solver schemes; dimensionless profiles |
| `src/types.jl` | `SurfaceFluxConfig`, `FluxSpecs`, `SolverOptions`, `SurfaceFluxConditions`, moisture/gustiness model types |
| `src/bulk_fluxes.jl` | Sensible/latent heat, evaporation, buoyancy, momentum fluxes, bulk Richardson number |
| `src/exchange_coefficients.jl` | Drag/heat exchange coefficients and conductance |
| `src/physical_scales.jl` | u\*, θ\*, q\*, variances, Obukhov length and stability parameter |
| `src/roughness_lengths.jl` | Constant, COARE 3.0, Raupach roughness models; combined u\*–roughness solver |
| `src/wind_and_gustiness.jl` | Effective wind speed and gustiness (constant, Deardorff) |
| `src/profile_recovery.jl` | `compute_profile_value` for diagnosing variables at arbitrary heights |
| `src/utilities.jl` | `surface_density`, geopotential helpers, `non_zero`, Gauss-Legendre quadrature |
| `src/input_builders.jl` | `build_surface_flux_inputs`: normalizes user inputs into a NamedTuple |
| `src/Parameters.jl` | `SurfaceFluxesParameters` and accessors |
| `ext/CreateParametersExt.jl` | ClimaParams-based constructors (weak dependency) |
| `test/` | Test suite (`runtests.jl`; GPU via `runtests_gpu.jl`) |
| `docs/` | Documentation source (`docs/src/`) and shared dev-guides (`docs/dev-guides/`) |

## Local norms

- For package tests, prefer `Pkg.test()` over manually `include`ing `test/runtests.jl`, so test-only dependencies load through the package test path.
- Physics code is dimensional: carry SI units in docstrings (square brackets, e.g. `[W/m^2]`, `[kg/kg]`, `[m/s]`). Keep sign conventions explicit (fluxes are positive upward).
- Match existing style: explicit names, narrow imports, comments that explain *why*. Unicode variable names (`ζ`, `θ`, `ϕ`, `ψ`, `Ψ`, `κ`, `ρ`) match the math.
- Docstrings follow [docs/dev-guides/code-quality/documentation_policy.md](docs/dev-guides/code-quality/documentation_policy.md): an indented signature line, single-`#` section headings in plural standard form (`# Arguments`, `# Returns`, `# Fields`, `# Examples`, `# Notes`), and `[`name`](@ref)` cross-references for every type/function mentioned. Use `raw"""..."""` for docstrings with LaTeX backslashes.
- Run `julia -e 'using JuliaFormatter; format(".")'` before committing code (config in `.JuliaFormatter.toml`, margin 92).

## Self-correction

- If the source layout table above is discovered to be stale, update it.
- If the user gives a correction about how work should be done in this repo, add it to `Local norms` or another clearly labeled persistent section in this file so future sessions inherit it.
