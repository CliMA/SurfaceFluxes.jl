[v1.0.1] Documentation audit and cleanup. Export `compute_profile_value` (previously
documented but not exported). Fix docstring/Documenter rendering issues (a stray tab in the
`windspeed` math block, and `[units] (...)` patterns that Documenter misparsed as links;
resolve the `build_surface_flux_inputs` cross-references). Generalize the variance universal
functions (`MomentumVariance`/`HeatVariance`) to the abstract parameter type, removing
duplicated per-parameterization placeholders. Make `compute_theta_star`/`compute_q_star`
fall back to interior values when the surface guess is `nothing` (previously they could
error). Vendor the shared CliMA DeveloperGuides under `docs/dev-guides/` with a monthly
auto-sync workflow, and add `AGENTS.md`.

**Behavior changes:**

- The default roughness used by `surface_fluxes` when `config` is omitted is now
  `z0m = 2e-4` m, `z0s = 2e-5` m (the `ConstantRoughnessParams` keyword defaults), unified
  from the previous `1e-3`/`1e-3`. Callers that pass an explicit `config` (including all
  regression tests) are unaffected.
- `u_variance` is **renamed to `surface_tke`** (it returns the surface-layer turbulent kinetic
  energy, not the streamwise variance `σ_u²`); `u_variance` is retained as a deprecated alias.
  It also now uses the TKE-based similarity
  (Tan et al. 2018) for the `GryanikParams` and `GrachevParams` parameterizations as well;
  previously those silently fell back to the streamwise-variance (Panofsky et al. 1977) form,
  ignoring the convective velocity scale. `BusingerParams` results are unchanged.

[PR 230/231] Updates docs; minor bug fix; additional tests. Release of v1.0

[PR 212] Refactor of SurfaceFluxes.jl: Consistently use stability parameter in all solvers and as inputs to many functions. Added functionality for wind speed dependent roughness lengths (Charnock, COARE3) and option to use functions to compute surface temperature/humidity. 

[PR 211] Fixes a bug in sensible heat flux

[PR 206] Updates UniversalFunctions.jl: Update functions to avoid catastrophic cancellations. Add continuity and linearisation tests + fix bug in the near-neutral limit. 

[PR 186] Removes unused stability function types (Holtslag, Cheng, Beljaars). Currently supports (Businger, Grachev, Gryanik) types.  
