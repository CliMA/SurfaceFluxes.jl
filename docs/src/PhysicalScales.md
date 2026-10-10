# Physical Scales

In Monin-Obukhov Similarity Theory (MOST), turbulent fluxes are parameterized using characteristic physical scales. These scales represent the turbulent fluctuations of velocity and scalars in the surface layer.

## Friction Velocity ($u_*$)

The friction velocity $u_*$ is the characteristic velocity scale of the turbulence, related to the surface kinematic momentum flux (stress) $\tau/\rho$:

```math
 u_*^2 = \frac{|\tau|}{\rho} = \left( (\overline{u'w'})^2 + (\overline{v'w'})^2 \right)^{1/2}.
```

In `SurfaceFluxes.jl`, this is computed by [`compute_ustar`](@ref). $u_*$ is derived from the effective wind speed difference $\Delta U = U_{\text{eff}}$ (gustiness included) using the dimensionless momentum profile $F_m$,

```math
 u_* = \frac{\kappa \Delta U}{F_m(\zeta, ...)},
```

where $F_m = \ln(\Delta z_{\text{eff}}/z_{0m}) - \psi_m(\zeta) + \psi_m(\zeta z_{0m}/\Delta z_{\text{eff}})$ is the dimensionless profile integral for momentum (plus the roughness-sublayer correction, if configured). With gustiness ($U_{\text{eff}} > |\Delta\mathbf{u}|$), the stress the code returns is smaller than $\rho u_*^2$: $|\tau|/\rho = u_*^2 |\Delta\mathbf{u}| / U_{\text{eff}}$.

The Obukhov length and the stability parameter follow from $u_*$ and the surface buoyancy flux $B$: $L = -u_*^3/(\kappa B)$ ([`obukhov_length`](@ref)) and $\zeta = (\Delta z - d)/L$.

## Scalar Scales ($\theta_*, q_*$)

Similar scales are defined for potential temperature ($\theta$) and specific humidity ($q$).

**Temperature Scale ($\theta_*$):**
Related to the kinematic potential temperature flux $\overline{w'\theta'}$:

```math
 u_* \theta_* = -\overline{w'\theta'}.
```

Computed as:

```math
 \theta_* = \frac{\kappa \Delta \theta}{F_h(\zeta, ...)},
```

with the interior-minus-surface difference $\Delta\theta = [c_{pd}(T_{\text{int}} - T_{\text{sfc}}) + g(\Delta z - d)]/c_{pm}$ of the dry static energy, divided by the moist heat capacity.

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

- **Neutral/Stable ($\zeta \ge 0$):** $\text{TKE} = 3.75 u_*^2$.
- **Unstable ($\zeta < 0$):** $\text{TKE} = 3.75 u_*^2 + 0.2 w_*^2 + u_*^2 (-\zeta)^{2/3}$, where the convective velocity scale $w_*$ depends on the boundary layer height $z_i$ (taken to be a fixed parameter) and the buoyancy flux implied by $\zeta$.

The streamwise variance $\sigma_u^2 = (u_* \phi_{\sigma u})^2$ (Panofsky et al. 1977) is available separately through the universal function `phi(uf, ζ, MomentumVariance())`. Panofsky et al. express $\phi_{\sigma u}$ in terms of the mixed-layer stability $z_i/L$; this function is evaluated at the local $\zeta$.

### Scalar Variance ($\sigma_\phi^2$)

The variance of scalars (temperature, humidity), computed by [`scalar_variance`](@ref):

```julia
scalar_variance(param_set, scale, ζ)
```

With the scalar scale $s_*$ (for example, $\theta_*$ or $q_*$):

- **Stable:** $\sigma_s^2 = 4 s_*^2$.
- **Unstable:** $\sigma_s^2 = 4 s_*^2 (1 - 8.3\zeta)^{-2/3}$ (Tan et al. 2018), whose free-convection limit $\sigma_s \approx 0.99 |s_*| (-\zeta)^{-1/3}$ matches [Wyngaard et al. (1971)](https://doi.org/10.1175/1520-0469(1971)028<1171:LFCSAT>2.0.CO;2).

## References

- Panofsky, H. A., Tennekes, H., Lenschow, D. H., & Wyngaard, J. C. (1977). The characteristics of turbulent velocity components in the surface layer under convective conditions. *Boundary-Layer Meteorology*, 11, 355–361. [DOI: 10.1007/BF02186086](https://doi.org/10.1007/BF02186086)
- Tan, Z., Kaul, C. M., Pressel, K. G., Cohen, Y., Schneider, T., & Teixeira, J. (2018). An extended eddy-diffusivity mass-flux scheme for unified representation of subgrid-scale turbulence and convection. *Journal of Advances in Modeling Earth Systems*, 10, 770–800. [DOI: 10.1002/2017MS001162](https://doi.org/10.1002/2017MS001162)
- Wyngaard, J. C., Coté, O. R., & Izumi, Y. (1971). Local free convection, similarity, and the budgets of shear stress and heat flux. *Journal of the Atmospheric Sciences*, 28, 1171–1182. [DOI: 10.1175/1520-0469(1971)028<1171:LFCSAT>2.0.CO;2](https://doi.org/10.1175/1520-0469(1971)028<1171:LFCSAT>2.0.CO;2)
