[v1.5.0] The surface state applies at the displacement height: the geopotential of the
surface temperature and humidity is `Φ_sfc + g d`
(`surface_geopotential(param_set, inputs)`), so the dry static energy difference that
drives the sensible heat flux is `cp (T_int - T_sfc) + g (Δz - d)`, and the hydrostatic
extrapolation of the surface density spans `Δz - d`
(`surface_density(param_set, inputs, T_sfc, q_vap_sfc)`). Fluxes over a canopy no longer
include the air column below the displacement height: results change for every `d ≠ 0`
(by `g d / cp` in the temperature difference, 0.15 K for `d = 15` m), are unchanged for
`d = 0`, and are unchanged when the displacement height and the reference level are
raised together. `surface_geopotential` takes the parameter set as its first argument;
the one-argument form is deprecated and returns the geopotential of the ground.

[v1.5.0] Reference level conventions: `SurfaceFluxConfig` has the new field
`reference_level`, `ReferenceAboveSurface()` by default (`Δz` measured from the surface)
or `ReferenceAboveApparentSink()` (`Δz` measured from `d + z0m`, as in the Community Land
Model). The solver and `screen_level_values` convert the latter to the former with
`reference_above_surface(param_set, inputs)`, so a forcing height below a tall canopy
remains valid. The conversion requires a roughness model independent of `u★`.

[v1.5.0] Raupach (1994) canopy roughness: `displacement_height(spec, roughness_inputs)` returns
the zero-plane displacement height `d` of the canopy from the same inputs as the roughness
length, for callers to pass as the input `d` of the solve and to place sub-canopy and
screen-level diagnostics. The roughness input `PAI` is the plant area index `Λ` (leaves plus
stems, `LAI + SAI`; the field name `LAI` is a deprecated alias, to be removed in the next
breaking release), and the drag partition uses the frontal area index `λ`, the area facing the wind
per unit ground area. `RaupachRoughnessParams` has the new fields `frontal_area_ratio`
(`λ / Λ`, 0.5 for isotropically oriented elements; previously fixed), `λ_min` (floor on `λ`,
zero by default; a positive floor stands in for stems and branches when the input counts
leaves only), and the constants `ustar_Uh_max` (0.3) and `c_w` (2, the roughness-sublayer
depth ratio, which sets the influence function `Ψ_h = ln c_w - 1 + 1 / c_w = 0.193`).
`RaupachRoughnessParams(toml_dict)` reads all coefficients from ClimaParams (v1.3.1 or
later), as the `raupach_*` parameters. The momentum
roughness length is bounded below by `z0m_fixed`. `frontal_area_ratio` enters the drag partition (Eq. 7) only; the displacement
height (Eq. 8) is an empirical fit in the plant area index. The scalar helpers
`frontal_area_index`, `canopy_area_index`, `raupach_displacement_fraction`, and
`raupach_roughness_fraction` expose the closed forms;
an excess resistance `kB⁻¹` is set through `stanton_number = exp(-kB⁻¹)`.

[v1.5.0] Raupach (1994) correction: the displacement height (Eq. 8) depends on the canopy
area index `Λ = 2λ`, not on the frontal area index `λ` as before. `d / h` is larger (0.66
instead of 0.56 at `Λ = 1`, matching Fig. 1b of the paper), and `z0m / h`, proportional to
`1 - d / h`, is smaller (peak 0.091 instead of 0.114).

[v1.5.0] Screen-level reconstruction `screen_level_values(param_set, sc, inputs, z_screen,
z_anemometer)`: the air temperature and vapor specific humidity at a height above the
apparent sink for heat `d + z0h`, and the wind speed at a height above the apparent sink
for momentum `d + z0m`, from the Monin-Obukhov profiles of a solve. Scalars follow
`X_sfc + (X_int - X_sfc) F̂_h(z) / F̂_h(Δz_eff)` at the effective Obukhov length `L_eff`, so
the profile reproduces the fluxes also under a stability cap; the temperature follows the
dry static energy, which the sensible heat flux is computed from; the wind speed is
`u* F̂_m(z) / κ`, the effective wind speed of the solve at the reference level. Heights are
clamped between the roughness length and the reference level (`dimensionless_profile_value`).
The screen and anemometer values are point values of the profiles; after a layer-average
solve, the reference profile follows `LayerAverageScheme`.

[v1.5.0] `SurfaceFluxConditions` has the new field `ζ_eff = min(ζ, ζ_cap) = Δz_eff / L_eff`,
the stability parameter at which the exchange coefficients were evaluated. The positional
constructors without it derive it from `ζ`, `L_MO`, and `L_eff`.

[v1.5.0] Reference level validity: `surface_fluxes` returns `NaN` in every field with
`converged = false` when the reference level lies at or below a roughness length
(`Δz - d <= max(z0m, z0h)`, `reference_height_valid`), in every solver mode and without
throwing, so that it runs in kernels; `check_reference_height(Δz, d, z0m, z0h)` throws an
`ArgumentError` for host-side validation of a configuration.

[v1.5.0] Gustiness floor accessors `minimum_wind_speed(spec, param_set)`, the minimum effective wind
speed a gustiness model imposes (zero for models without a floor), and
`without_floor(spec)`, the same model with a zero floor, for callers that fold the floor
into the wind they pass to the solve.

[v1.5.0] New helper `interior_vapor_specific_humidity(inputs)`, the interior total specific
humidity less the condensate, used by the evaporation and bulk Richardson number.

[v1.4.0] New gustiness model `FlooredDeardorffGustinessSpec(u_min)`: the larger of a minimum
wind speed and the Deardorff convective gustiness `β w*`, with the convective part
evaluated in closed form within the stability solve. At a stability parameter `ζ`, the
bulk relations `u* = κ U / F_m(ζ)` and `θv* = κ Δθv / F_h(ζ)` make the buoyancy flux
linear in the effective wind speed `U`, so `U = β w*(U)` has the solution
`U² = β³ κ² (g/θv) z_i Δθv / (F_m F_h)` (`free_convection_wind_speed`), with the profile
integrals of the solver's discretization scheme. The gustiness
does not depend on `u*` (`depends_on_ustar` is `false`), so the friction velocity follows
from `ζ` without the inner Brent solve, and the free-convection limit is well posed at
every `ζ`; `DeardorffGustinessSpec` evaluates the gustiness from the buoyancy flux
implied by `ζ` and the current `u*`, which has no consistent `u*` beyond the
free-convection limit. New helper `virtual_pottemps` returns the surface and interior
virtual potential temperatures used by `state_bulk_richardson_number` and the closed
form.

[v1.3.0] Raupach (1994) momentum roughness: `u★ / U(h)` is capped at 0.3, the sheltering
limit of Eq. 7. `z0m / h` now peaks at `λ ≈ 0.29` (`LAI ≈ 0.58`, `z0m / h ≈ 0.11`) and
decreases for denser canopies; the uncapped form kept increasing with `LAI`.

[v1.3.0] Friction velocity solve with `ustar`-dependent gustiness or roughness:

- When the bracket `ustar ∈ [1e-4, 4]` m/s of the inner Brent solve contains no
  consistent `ustar`, the endpoint on the side of the root is returned (previously the
  endpoint of smaller residual, `1e-4` m/s). With `DeardorffGustinessSpec`, the gustiness
  at fixed `ζ` is proportional to `ustar`, and no consistent `ustar` exists for `ζ` more
  unstable than the free-convection limit. Returning `1e-4` m/s there let the ζ solve
  converge to spurious roots with a vanishing friction velocity and a large sensible heat
  flux over rough surfaces (`z0m` of order 1 m) in unstable conditions.

[v1.3.0] Stability caps for stable stratification:

- New `SurfaceFluxConfig` field `stability_cap` (fifth positional argument; default
  `NoStabilityCap()`, standard MOST). `ConstantStabilityCap(ζ_max)` caps the stability
  parameter in the flux-profile relations at a constant;
  `MaxHeatFluxStabilityCap()` caps it at the stability `ζ_p` of maximum sensible heat flux
  at fixed wind speed (`max_heat_flux_stability`), beyond which MOST predicts a heat flux
  that decreases with increasing stratification (runaway cooling and decoupling).
- Beyond the cap, the exchange coefficients and similarity scales are held at their
  values at the cap, so the bulk Richardson number is linear in `ζ` and the MOST solve
  has a root (for caps within the solver's range `ζ ≤ 100`, roots beyond `|ζ| = 100`
  are bracketed by an extended probe). The
  returned `ζ` and `L_MO` are those implied by the fluxes. The cap also applies to the
  diagnostic heat conductance when fluxes are prescribed, and to the conductance seen by
  surface-state callbacks through `heat_conductance`.
- `SurfaceFluxConditions` has a new field `L_eff` (after `L_MO`), the effective Obukhov
  length `Δz_eff / min(ζ, ζ_cap)`. Pass it to `compute_profile_value` for profiles
  consistent with the capped fluxes. It equals `L_MO` without an active cap. The
  positional constructor accepts the fields with or without `L_eff` (then `L_eff = L_MO`).
- `ConstantStabilityCap` requires `ζ_max > 0`, is converted to the floating-point type
  of the parameter set, and is differentiable with respect to `ζ_max`. `heat_conductance`, `compute_ustar`, `compute_theta_star`, and `compute_q_star`
  compute the cap from `inputs.stability_cap` when `inputs.ζ_cap` is `nothing` (inputs from
  `build_surface_flux_inputs`), so they agree with `surface_fluxes` for the same
  configuration.

[v1.3.0] Roughness sublayer (RSL) corrections reworked (the RSL models have not yet been
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
