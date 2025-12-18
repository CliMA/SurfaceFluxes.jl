# Prescribed Conditions

SurfaceFluxes.jl supports flexible operating modes for specifying boundary conditions. While the default behavior is to solve for the surface fluxes iteratively using Monin-Obukhov Similarity Theory (MOST) given the temperature and humidity at the surface, users can also prescribe specific fluxes or exchange coefficients. It is also possible to compute the surface temperature and humidity interactively via callbacks. This is particularly useful for:
- Driving models with observed fluxes.
- Coupling to components that compute their own exchange coefficients (e.g., wave models).
- Determining skin temperature over land/ice and implementing moisture limitations on evaporation.
- Idealized experiments with fixed boundary conditions.

## Flux Specification

The operating mode is controlled by the optional `flux_specs` argument in [`surface_fluxes`](@ref). This argument accepts a [`FluxSpecs`](@ref) struct, which acts as a container for prescribed quantities:

```julia
FluxSpecs(shf=..., lhf=..., ustar=..., Cd=..., Ch=...)
```

| Field   | Symbol | Description | Units |
|:--------|:-------|:------------|:------|
| `shf`   | $H$    | Sensible Heat Flux | $\mathrm{W~m^{-2}}$ |
| `lhf`   | $LE$   | Latent Heat Flux | $\mathrm{W~m^{-2}}$ |
| `ustar` | $u_*$  | Friction Velocity | $\mathrm{m~s^{-1}}$ |
| `Cd`    | $C_d$  | Momentum Exchange Coefficient | - |
| `Ch`    | $C_h$  | Heat Exchange Coefficient | - |

The solver detects which combination of parameters is provided and dispatches to the appropriate routine.

## Operating Modes

### 1. Iterative Solver (Standard MOST)
- **Inputs:** Surface state ($T_s, q_s, \mathbf{u}_s$) and Atmospheric state ($T_a, q_a, \mathbf{u}_a, z, d$).
- **Unknowns:** Fluxes ($H, LE, \boldsymbol{\tau}$), Coefficients ($C_d, C_h$), Stability ($\zeta$).

This is the default mode when `flux_specs` is empty or `nothing`. The solver iterates to find the Monin-Obukhov stability parameter $\zeta$ that satisfies the similarity relations.

```julia
# Standard call (positional arguments)
conditions = surface_fluxes(param_set, T_int, ..., u_sfc, ...)
```

!!! note "Simplified Examples"
    The examples below use simplified syntax (e.g., `flux_specs=...`) for clarity. In the actual API, [`surface_fluxes`](@ref) uses positional arguments. Users must provide all preceding arguments or use `nothing` for optional inputs. See the [API Reference](API.md) for the exact signature.

### 2. Prescribed Coefficients
- **Inputs:** Coefficients ($C_d, C_h$).
- **Unknowns:** Fluxes ($H, LE, \boldsymbol{\tau}$).

When both `Cd` and `Ch` are prescribed, the fluxes are computed directly from the bulk aerodynamic formulas without iteration. This is typical when coefficients are diagnosed by an external parameterization (e.g., a wave model).

```julia
specs = FluxSpecs(Cd = 1.2e-3, Ch = 1.2e-3)
conditions = surface_fluxes(param_set, ..., flux_specs=specs)
```

### 3. Fully Prescribed Fluxes
- **Inputs:** Fluxes ($H, LE$) and Friction Velocity ($u_*$).
- **Unknowns:** Stability ($\zeta$), Coefficients ($C_d, C_h$).

When `shf`, `lhf`, and `ustar` are all provided, the solver bypasses the flux calculation. It uses the prescribed values to:
1.  Compute the Monin-Obukhov length $L$.
2.  Diagnose the stability parameter $\zeta = z/L$.
3.  Back-calculate the exchange coefficients consistent with these fluxes.

```julia
specs = FluxSpecs(shf = 20.0, lhf = 100.0, ustar = 0.3)
conditions = surface_fluxes(param_set, ..., flux_specs=specs)
```

### 4. Prescribed Heat Fluxes and Drag Coefficient
- **Inputs:** Fluxes ($H, LE$) and Drag Coefficient ($C_d$).
- **Unknowns:** Friction Velocity ($u_*$), Stability ($\zeta$).

When `shf`, `lhf`, and `Cd` are provided (but not `ustar`):
1.  Friction velocity is derived from the drag law: $u_* = \sqrt{C_d} U_{\text{eff}}$.
2.  Stability $\zeta$ is computed from the fluxes and the derived $u_*$.

```julia
specs = FluxSpecs(shf = 50.0, lhf = 150.0, Cd = 1.5e-3)
conditions = surface_fluxes(param_set, ..., flux_specs=specs)
```

## State Callbacks

In coupled systems, the surface state (e.g., skin temperature $T_s$) may respond rapidly to the surface fluxes. To ensure consistency without external fixed-point iteration, [`surface_fluxes`](@ref) accepts callback functions that update the surface state *within* the MOST iteration loop.

```julia
# Update Ts based on current flux guess
function update_Ts(T_sfc, local_geometry, surface_fluxes_conditions)
    # ... physics to update skin temperature ...
    return new_T_sfc
end

conditions = surface_fluxes(..., update_T_sfc=update_Ts)
```

Available callbacks:
- `update_T_sfc`: Updates surface temperature.
- `update_q_vap_sfc`: Updates surface specific humidity.

This mechanism ensures that the final fluxes and surface state are in equilibrium with respect to the surface energy/moisture balance.