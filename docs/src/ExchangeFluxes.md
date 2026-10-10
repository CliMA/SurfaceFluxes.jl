# Exchange Fluxes

SurfaceFluxes.jl provides a robust interface for computing surface-atmosphere exchange. The fluxes are calculated using bulk aerodynamic formulas, parameterized by non-dimensional exchange coefficients derived from Monin-Obukhov Similarity Theory (MOST).

## Bulk Fluxes

The `SurfaceFluxes` module exports functions to compute the surface fluxes directly. These are the primary physical quantities of interest for coupling.

### 1. Momentum Fluxes ($\tau$)

The surface stress vector [$\mathrm{N~m^{-2}}$], representing the transfer of momentum from the atmosphere to the surface.

```julia
momentum_fluxes(Cd, inputs, ρ_sfc, gustiness)
```

**Formula:**

```math
\tau_{x,y} = -\rho_{\text{sfc}} C_d U_{\text{eff}} \Delta u_{x,y}
```

where:

- The drag coefficient is denoted by $C_d$.
- The effective wind speed $U_{\text{eff}}$ includes gustiness effects.
- The wind speed component differences are $\Delta u_{x,y}$.
- The surface density $\rho_{\text{sfc}}$ is computed internally by extrapolating the interior
  pressure hydrostatically over $\Delta z - d$ to the displacement height, where the surface
  state applies.

### 2. Evaporation ($E$)

The mass flux of water vapor [$\mathrm{kg/m^2/s}$], driven by the specific humidity gradient.

```julia
evaporation(param_set, inputs, g_h, q_vap_int, q_vap_sfc, ρ_sfc, model)
```

**Formula:**

```math
E = -\rho_{\text{sfc}} g_h (q_{\text{vap,int}} - q_{\text{vap,sfc}})
```

where:

- The surface air density is $\rho_{\text{sfc}}$.
- The term $g_h$ represents the **conductance** for heat/scalars [$\mathrm{m~s^{-1}}$].
- The term $(q_{\text{vap,int}} - q_{\text{vap,sfc}})$ is the specific humidity difference between the interior (atmosphere) and the surface.

If a latent heat flux is prescribed in `inputs`, the evaporation is derived from it: $E = \text{LHF} / L_{v,0}$.

### 3. Latent Heat Flux (LHF)

The energy flux associated with the phase change of water [$\mathrm{W/m^2}$].

```julia
latent_heat_flux(param_set, inputs, E, model)
```

**Formula:**

```math
\text{LHF} = L_{v,0} E
```

where $L_{v,0}$ is the latent heat of vaporization at the reference temperature.

### 4. Sensible Heat Flux (SHF)

The energy flux driven by the temperature difference [$\mathrm{W/m^2}$].

```julia
sensible_heat_flux(param_set, inputs, g_h, T_int, T_sfc, ρ_sfc, E)
```

**Formula:**

```math
\text{SHF} = -\rho_{\text{sfc}} g_h (\text{DSE}_{\text{int}} - \text{DSE}_{\text{sfc}}) + \text{VSE}_{\text{sfc}} \times E
```

This formulation accounts for the enthalpy transport due to sensible heat transfer:

1. **Dry Static Energy Term**: $-\rho_{\text{sfc}} g_h \Delta \text{DSE}$. Driven by the dry static energy (potential temperature) gradient. $\text{DSE}_{\text{sfc}}$ is evaluated at $T_{\text{sfc}}$ and the geopotential $\Phi_{\text{sfc}} + g d$ of the surface state, so $\Delta\text{DSE} = c_{pd}(T_{\text{int}} - T_{\text{sfc}}) + g(\Delta z - d)$ (see [Reference Level](SurfaceFluxes.md#Reference-Level)).
2. **Mass Transfer Term**: $\text{VSE}_{\text{sfc}} \times E$. The vapor static energy $\text{VSE}_{\text{sfc}} = c_{pv}(T_{\text{sfc}} - T_0) + \Phi_{\text{sfc}} + g d$ is the sensible enthalpy plus potential energy carried by the evaporated water; its latent heat is in LHF.

See [Yatunin et al. (2026)](https://doi.org/10.1029/2025MS005014) for a derivation of these formulas and a detailed explanation of how they result in an energetically consistent formulation.

## Exchange Coefficients

The non-dimensional exchange coefficients relate the fluxes to the bulk gradients. They are derived from the integrated similarity profiles ($F_m, F_h$) computed by the `UniversalFunctions` module.

### Drag Coefficient ($C_d$)

For momentum exchange:

```math
C_d = \left( \frac{\kappa}{F_m(\Delta z_{\text{eff}}, \zeta, z_{0m})} \right)^2
```

Momentum flux: $\boldsymbol\tau = -\rho_{\text{sfc}} C_d U_{\text{eff}} \Delta\mathbf{u}$ (see above). Here and in $C_h$, $F$ is the RSL-corrected profile $\widehat F = F + P$ when a roughness-sublayer model is configured, evaluated at the capped $\zeta$ when a stability cap is set.

### Heat Exchange Coefficient ($C_h$)

For heat and scalar exchange:

```math
C_h = \frac{\kappa^2}{F_m(\Delta z_{\text{eff}}, \zeta, z_{0m}) F_h(\Delta z_{\text{eff}}, \zeta, z_{0h})}
```

### Conductance ($g_h$)

The bulk scalar fluxes are proportional to the **conductance** $g_h$, which is the product of the non-dimensional heat exchange coefficient $C_h$ and the effective wind speed:

```math
g_h = C_h U_{\text{eff}}
```

## Roughness Length Models

The surface roughness lengths ($z_{0m}, z_{0h}$) parameterize the effect of surface irregularities on the wind and scalar profiles. SurfaceFluxes.jl supports several models via [`SurfaceFluxConfig`](@ref).

### Constant Roughness

`ConstantRoughnessParams` uses fixed values for $z_{0m}$ and $z_{0h}$. This is typical for static land surfaces where roughness is prescribed.

### COARE 3.0 (Ocean)

`COARE3RoughnessParams` implements the COARE 3.0 algorithm ([Fairall et al., 2003](https://doi.org/10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2)) for open ocean surfaces.

- **Momentum Roughness ($z_{0m}$)**: A sum of a smooth flow limit (viscous) and a rough flow limit (Charnock relationship, with a Charnock parameter $\alpha$ that varies with wind speed):

```math
z_{0m} = 0.11 \frac{\nu}{u_*} + \alpha \frac{u_*^2}{g}
```

- **Scalar Roughness ($z_{0h}$)**: Based on Reynolds number scaling, $z_{0h} = \min(1.1\times10^{-4}\ \mathrm{m},\ 5.5\times10^{-5}\ \mathrm{m}\ Re_*^{-0.6})$ with $Re_* = z_{0m} u_* / \nu$ (Fairall et al. 2003).

### Raupach (Land/Canopy)

`RaupachRoughnessParams` implements the [Raupach (1994)](https://doi.org/10.1007/BF00709229) model for vegetation canopies.

- Estimates $z_{0m}$ and the displacement height $d$ from the canopy height $h$ and the plant area index `PAI`, the one-sided area of leaves, stems, and branches per unit ground area: the sum of the leaf and stem area indices, LAI + SAI. (The field name `LAI` is accepted as a deprecated alias for `PAI`.) The displacement height depends on `PAI` directly (raised to $\lambda_{min}/f_\lambda$ where the floor below is active). The split of the drag between the ground and the plants depends on the frontal area index $\lambda = \max(f_\lambda\,\mathrm{PAI}, \lambda_{min})$, the area the canopy presents to the wind per unit ground area. The frontal area ratio $f_\lambda$ (`frontal_area_ratio`) is 0.5 by default, the value for randomly oriented elements. The optional floor $\lambda_{min}$ (zero by default) keeps a canopy rough when the input counts leaves only, as for a leafless deciduous forest. $z_{0m}$ is at least the fixed roughness length `z0m_fixed`.
- [`surface_fluxes`](@ref) takes $d$ as a separate argument. [`SurfaceFluxes.displacement_height`](@ref) computes it from the same canopy inputs as $z_{0m}$, so that the two are consistent.
- The scalar roughness length is $z_{0h} = r\,z_{0m}$, with the fixed ratio $r = \exp(-kB^{-1})$ (field `stanton_number`, 0.1 by default, so $kB^{-1} \approx 2.3$).
- Useful for dynamic vegetation models.

## Gustiness

In unstable conditions, especially when the mean wind speed approaches zero (free convection limit), convective eddies generated by surface heating maintain turbulent exchange. Standard bulk formulas using only the mean wind speed difference ($\Delta U$) would erroneously predict zero fluxes.

To account for this, `SurfaceFluxes.jl` uses an **effective wind speed** ($U_{\text{eff}}$):

```math
U_{\text{eff}} = \max\left( \left(\Delta u^2 + \Delta v^2\right)^{1/2}, \; U_{\text{gust}} \right)
```

where $U_{\text{gust}}$ is a parameterized gustiness velocity scale representing the contribution of sub-grid eddies.

### Parameterizations and Dispatch

The gustiness formulation is controlled by the `gustiness` field of [`SurfaceFluxConfig`](@ref), a subtype of `AbstractGustinessSpec`, and [`surface_fluxes`](@ref) dispatches on its type. Land and ocean points can use different configurations; within one broadcast, a single concrete specification type keeps the computation type-stable.

#### 1. Constant Gustiness

```julia
ConstantGustinessSpec(value)
```

Uses a fixed tuning parameter, e.g., $U_{\text{gust}} = 1.0 \, \mathrm{m~s^{-1}}$. Often used over land surfaces.

#### 2. Deardorff Gustiness

```julia
DeardorffGustinessSpec()
```

Scales $U_{\text{gust}}$ with the convective velocity scale $w_*$, following [Deardorff (1970)](https://doi.org/10.1175/1520-0469(1970)027<1211:CVATSF>2.0.CO;2) and [Beljaars (1995)](https://doi.org/10.1002/qj.49712152203). Beljaars (1995) adds $\beta w_*$ to the mean wind in quadrature; SurfaceFluxes.jl takes the larger of the two. This is physically robust for the unstable boundary layer over fluid surfaces (ocean/lakes).

```math
U_{\text{gust}} = \beta w_* = \beta (B z_i)^{1/3}
```

where:

- The surface buoyancy flux is denoted by $B$.
- The boundary layer height is $z_i$, the fixed parameter `gustiness_zi` (1000 m by default).
- The scaling coefficient is $\beta$ (`gustiness_coeff`, 1.25 by default).

Since $B$ depends on the fluxes, and the fluxes depend on $U_{\text{eff}}$ (and thus $B$), this introduces a nonlinear coupling that is resolved by an iterative solver. Within the stability solve, $B$ is evaluated from the current friction velocity, $B = -u_*^3 \zeta / [\kappa (\Delta z - d)]$, so that at fixed $\zeta$ the gustiness is proportional to $u_*$. Beyond the free-convection limit, no $u_*$ is consistent with this gustiness, and the solve settles where one exists.

#### 3. Floored Deardorff Gustiness

```julia
FlooredDeardorffGustinessSpec(u_min)
```

The larger of a minimum wind speed $u_{\min}$ and the Deardorff gustiness, with the convective part evaluated in closed form from the surface and atmospheric state. At a stability parameter $\zeta$, the bulk relations $u_* = \kappa U_{\text{eff}} / F_m(\zeta)$ and $\theta_{v*} = \kappa \Delta\theta_v / F_h(\zeta)$ (defined with the surface excess, so that $u_* \theta_{v*} = \overline{w'\theta_v'}$) make the buoyancy flux $B = (g/\theta_v) u_* \theta_{v*}$ linear in $U_{\text{eff}}$, so that $U_{\text{eff}} = \beta w_*(U_{\text{eff}})$ has the solution

```math
U_{\text{eff}}^2 = \beta^3 \kappa^2 \frac{g}{\theta_v} z_i \frac{\Delta\theta_v}{F_m(\zeta) F_h(\zeta)},
```

where $\Delta\theta_v$ is the virtual potential temperature excess of the surface over the air, $\theta_v$ is the air value, and $F_m$, $F_h$ are the dimensionless profile integrals of momentum and heat. The convective part vanishes when the surface is not warmer than the air, where the floor applies. Because this gustiness does not depend on $u_*$, with a roughness model independent of $u_*$ the friction velocity follows from $\zeta$ in closed form, and the free-convection limit is well posed at every $\zeta$. `FlooredDeardorffGustinessSpec(0)` is the pure convective gustiness in this form. This model is suited to land surfaces, where the surface temperature responds quickly to the fluxes and a small floor keeps the exchange from vanishing in calm, stable conditions.

#### Minimum Wind Speed

[`minimum_wind_speed`](@ref)`(spec, param_set)` returns the lowest effective wind speed a gustiness model allows: the value of a `ConstantGustinessSpec`, the floor $u_{\min}$ of a `FlooredDeardorffGustinessSpec`, and zero for a `DeardorffGustinessSpec`, whose gustiness vanishes in stable conditions. [`without_floor`](@ref)`(spec)` returns the same model with the floor set to zero. A model that applies the floor to the wind before passing it to the solve, for example a canopy model that reduces the wind above the canopy to the wind at the ground, passes `without_floor(spec)` so that the floor is not applied twice.

## Reference

- Beljaars, A. C. M. (1995). The parametrization of surface fluxes in large-scale models under free convection. *Quarterly Journal of the Royal Meteorological Society*, 121, 255-270. [DOI: 10.1002/qj.49712152203](https://doi.org/10.1002/qj.49712152203)

- Deardorff, J. W. (1970). Convective velocity and temperature scales for the unstable planetary boundary layer and for Rayleigh convection. *Journal of the Atmospheric Sciences*, 27, 1211-1213. [DOI: 10.1175/1520-0469(1970)027<1211:CVATSF>2.0.CO;2](https://doi.org/10.1175/1520-0469(1970)027<1211:CVATSF>2.0.CO;2)

- Fairall, C. W., Bradley, E. F., Hare, J. E., Grachev, A. A., & Edson, J. B. (2003). Bulk parameterization of air–sea fluxes: Updates and verification for the COARE algorithm. *Journal of Climate*, 16, 571-591. [DOI: 10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2](https://doi.org/10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2)

- Raupach, M. R. (1994). Simplified expressions for vegetation roughness length and zero-plane displacement as functions of canopy height and area index. *Boundary-Layer Meteorology*, 71, 211-216. [DOI: 10.1007/BF00709229](https://doi.org/10.1007/BF00709229)

- Yatunin, D., Byrne, S., Kawczynski, C., Kandala, S., Bozzola, G., Sridhar, A., Shen, Z., Jaruga, A., Sloan, J., He, J., Huang, D.Z., Barra, V., Knoth, O., Ullrich, P., Schneider, T., 2026: The CliMA atmosphere dynamical core: Concepts, numerics, and scaling. *Journal of Advances in Modeling Earth Systems*, 18, e2025MS005014. [DOI:10.1029/2025MS005014](https://doi.org/10.1029/2025MS005014)
