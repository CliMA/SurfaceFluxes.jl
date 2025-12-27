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

### Velocity Variance ($\sigma_u^2$)
The variance of the horizontal wind speed components, computed by [`u_variance`](@ref):

```julia
u_variance(param_set, Δz_eff, ustar, ζ)
```

The parameterization depends on stability:
- **Neutral/Stable ($\zeta \ge 0$):** Proportional to $u_*^2$.
- **Unstable ($\zeta < 0$):** Includes a contribution from the convective velocity scale $w_*$, which depends on the boundary layer height $z_i$ (taken to be a fixed parameter) and the heat flux at the effective height $\Delta z_{\text{eff}}$.

### Scalar Variance ($\sigma_\phi^2$)
The variance of scalars (temperature, humidity), computed by [`scalar_variance`](@ref):

```julia
scalar_variance(param_set, scale, ζ)
```
- **Stable:** Proportional to $\phi_*^2$.
- **Unstable:** Scaling depends on $\zeta$ following standard similarity functions (e.g., [Wyngaard et al., 1971](https://doi.org/10.1175/1520-0469(1971)028<1171:LFCSAT>2.0.CO;2)).

