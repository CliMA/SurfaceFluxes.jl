[main] Roughness sublayer (RSL) corrections reworked (the RSL models have not yet been
released):

- `PhysickGarrattRSL` and `HarmanFinniganRSL` are renamed `LinearRSL` and
  `ExponentialRSL` (fields `c_m`, `c_h`, `z_RSL`). The exponential RSL factor
  `exp(-c(1 - z/z_RSL))` is the form of Garratt (1980) and Physick & Garratt (1995).
- The RSL correction is now the integral of `φ(z/L) (1 - μ(z))/z`, consistent with
  `φ̂ = φ μ`: it depends on stability, includes the neutral Prandtl number for scalars, and
  is layer-averaged for `LayerAverageScheme`. The corrected profiles satisfy `F̂ ≥ F`, so
  exchange coefficients stay positive and finite over tall canopies and in strongly
  unstable conditions (previously `F̂` could become negative).
- The corrected profiles are now anchored at the RSL top: they coincide with MOST above
  the RSL, so `z0` and `d` are the apparent canopy values (e.g., `z0 ≈ 0.1h`,
  `d ≈ 0.67h`). Previously, the profiles were anchored at `z0`, which double-counted the
  RSL effect when combined with apparent roughness lengths.
- `compute_theta_star`, `compute_q_star`, and `compute_profile_value` now include the RSL
  correction. `rsl_profile_correction(uf_params, rsl_model, Δz_eff, ζ, z0, transport, scheme)`
  has a new signature; `rsl_corrected_profile` returns the corrected profile `F̂`.
- Parameters are validated (`0 ≤ c < 1` for `LinearRSL`, `c ≥ 0`, `z_RSL ≥ 0`), and
  floating-point parameters are converted to the type of the inputs.


- When `update_T_sfc` or `update_q_vap_sfc` callbacks are supplied, the solver is routed through a `solve_stability_param_cb` function. On each solver call the `inputs`
  argument received by the callback will have `inputs.T_sfc_guess` and
  `inputs.q_vap_sfc_guess` set to the **previous iteration's returned values**,
  not the original guesses passed by the caller. Callbacks that read these fields
  as a starting point for their physics should be aware that the values evolve
  across solver iterations.

[v1.2.0] AD compatibility tests now cover Enzyme (forward and reverse) via
DifferentiationInterface, in addition to ForwardDiff. Derivatives of sensible heat
flux with respect to surface temperature are checked against central finite differences
across stable, near-neutral, and unstable regimes.

**Behavior changes:**

- When `update_T_sfc` or `update_q_vap_sfc` callbacks are supplied, the solver now
  uses an iterative residual functor (`IterativeResidualFunction`) that advances
  the surface-state guess across ζ iterations. On each solver call the `inputs`
  argument received by the callback will have `inputs.T_sfc_guess` and
  `inputs.q_vap_sfc_guess` set to the **previous iteration's returned values**,
  not the original guesses passed by the caller. Callbacks that read these fields
  as a starting point for their physics should be aware that the values evolve
  across solver iterations.

[v1.1.0] More robust solve for the stability parameter ζ. The unbracketed secant
iteration (which could take near-singular steps to |ζ| ≫ 100 in very stable conditions), is replaced by a branchless, bracketed solve, with a fixed number of residual evaluations.

**Behavior changes:**

- For supercritical `Ri_b` (no root within `|ζ| <= 100`), ζ now deterministically
  saturates at the limit of the correct stability branch with `converged = false`,
  instead of returning a clamped stray secant iterate that could lie on the wrong branch.
- Converged roots may differ from the secant solver's at the level of the solver
  tolerances; regression-test results are unchanged within their tolerances.

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
