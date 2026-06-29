# Physical Scales

In Monin-Obukhov Similarity Theory (MOST), turbulent fluxes are parameterized using characteristic physical scales. These scales represent the turbulent fluctuations of velocity and scalars in the surface layer.

## Friction Velocity ($u_*$)

The friction velocity $u_*$ is the characteristic velocity scale of the turbulence, related to the surface kinematic momentum flux (stress) $\tau/\rho$:

```math
 u_*^2 = \frac{|\tau|}{\rho} = \left( (\overline{u'w'})^2 + (\overline{v'w'})^2 \right)^{1/2}.
```

In `SurfaceFluxes.jl`, this is computed by [`compute_ustar`](@ref). $u_*$ is derived from the wind speed difference $\Delta U$ using the dimensionless momentum profile $\phi_m$,

```math
 u_* = \frac{\kappa \Delta U}{F_m(\zeta, ...)},
```

where $F_m$ is the integrated stability correction function for momentum.

## Scalar Scales ($\theta_*, q_*$)

Similar scales are defined for potential temperature ($\theta$) and specific humidity ($q$).

**Temperature Scale ($\theta_*$):**
Related to the kinematic potential temperature flux $\overline{w'\theta'}$:

```math
 u_* \theta_* = -\overline{w'\theta'}.
```

Computed as:

```math
 \theta_* = \frac{\kappa \Delta \theta}{F_h(\zeta, ...)}
```

In `SurfaceFluxes.jl`, this scale is computed by [`compute_theta_star`](@ref).

**Humidity Scale ($q_*$):**
Related to the kinematic specific humidity flux $\overline{w'q'}$ (evaporation):

```math
 u_* q_* = -\overline{w'q'}.
```

Computed similarly to $\theta_*$ using the same heat stability function $F_h$. See [`compute_q_star`](@ref).

## Variances

`SurfaceFluxes.jl` also provides functions to estimate the variances of turbulent fluctuations, which are useful for higher-order closure models or statistical analysis.

!!! warning "Range of validity"
    The variance and TKE functions are **convective surface-layer closures and are independent of the chosen flux-profile parameterization**: Businger, Gryanik, and Grachev all return the same values, because [Grachev et al. (2007)](https://doi.org/10.1007/s10546-007-9177-6) and [Gryanik et al. (2020)](https://doi.org/10.1175/JAS-D-19-0255.1) do not define variance similarity functions. On the **stable** side they reduce to constant (neutral) values. This is a crude approximation — observed $\sigma_u/u_*$ *increases* with stability rather than staying constant, and Monin–Obukhov scaling of the horizontal-velocity variances breaks down in stable stratification. Use the stable-side variances as rough estimates only; they are *not* validated stable-boundary-layer similarity, even when a stable-boundary-layer parameterization (Gryanik, Grachev) is selected for the fluxes.

### Turbulent Kinetic Energy

The function [`surface_tke`](@ref) returns the turbulent kinetic energy (TKE) following Tan et al. (2018):

```julia
surface_tke(param_set, Δz_eff, ustar, ζ)
```

Returns the TKE, $(u_* \phi)^2$. The parameterization depends on stability:

- **Neutral/Stable ($\zeta \ge 0$):** $\text{TKE} = 3.75\, u_*^2$.
- **Unstable ($\zeta < 0$):** $\text{TKE} = 3.75\, u_*^2 + 0.2\, w_*^2 + u_*^2 (-\zeta)^{2/3}$, where the convective velocity scale $w_*$ depends on the boundary layer height $z_i$ (taken to be a fixed parameter) and the buoyancy flux implied by $\zeta$.

The streamwise variance $\sigma_u^2 = (u_* \phi_{\sigma u})^2$ (Panofsky et al. 1977) is available separately through the universal function `phi(uf, ζ, MomentumVariance())`.

### Scalar Variance ($\sigma_\phi^2$)

The variance of scalars (temperature, humidity), computed by [`scalar_variance`](@ref):

```julia
scalar_variance(param_set, scale, ζ)
```

- **Stable:** Proportional to $\phi_*^2$.
- **Unstable:** Scaling depends on $\zeta$ following standard similarity functions (e.g., [Wyngaard et al., 1971](https://doi.org/10.1175/1520-0469(1971)028<1171:LFCSAT>2.0.CO;2)).
