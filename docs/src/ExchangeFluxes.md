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
- The surface density $\rho_{\text{sfc}}$ is computed internally by hydrostatically extrapolating
  the interior pressure to the surface.

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

1. **Dry Static Energy Term**: $-\rho_{\text{sfc}} g_h \Delta \text{DSE}$. Driven by the dry static energy (potential temperature) gradient.
2. **Mass Transfer Term**: $\text{VSE}_{\text{sfc}} \times E$. Represents the "dry" enthalpy carried by the evaporating water vapor leaving the surface ($\text{VSE}$ is Vapor Static Energy, or the dry static energy carried by water vapor).

See [Yatunin et al. (2026)](https://doi.org/10.1029/2025MS005014) for a derivation of these formulas and a detailed explanation of how they result in an energetically consistent formulation.

## Exchange Coefficients

The non-dimensional exchange coefficients relate the fluxes to the bulk gradients. They are derived from the integrated similarity profiles ($F_m, F_h$) computed by the `UniversalFunctions` module.

### Drag Coefficient ($C_d$)

For momentum exchange:

```math
C_d = \left( \frac{\kappa}{F_m(\Delta z_{\text{eff}}, \zeta, z_{0m})} \right)^2
```

Momentum flux: $\tau = \rho_{\text{sfc}} C_d U_{\text{eff}} \Delta U$.

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

The surface roughness lengths ($z_{0m}, z_{0h}$) parameterize the effect of surface irregularities on the wind and scalar profiles. SurfaceFluxes.jl supports several models via `SurfaceFluxFluxConfig`.

### Constant Roughness

`ConstantRoughnessParams` uses fixed values for $z_{0m}$ and $z_{0h}$. This is typical for static land surfaces where roughness is prescribed.

### COARE 3.0 (Ocean)

`COARE3RoughnessParams` implements the COARE 3.0 algorithm ([Fairall et al., 2003](https://doi.org/10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2)) for open ocean surfaces.

- **Momentum Roughness ($z_{0m}$)**: A sum of a smooth flow limit (viscous) and a rough flow limit (Charnock relationship, with a Charnock parameter $\alpha$ that varies with wind speed):

```math
z_{0m} = 0.11 \frac{\nu}{u_*} + \alpha \frac{u_*^2}{g}
```

- **Scalar Roughness ($z_{0h}$)**: Based on Reynolds number scaling.

### Raupach (Land/Canopy)

`RaupachRoughnessParams` implements the [Raupach (1994)](https://doi.org/10.1007/BF00709229) model for vegetation canopies.

- Calculates $z_{0m}$ and displacement height $d$ based on canopy height ($h$) and Leaf Area Index (LAI).
- Useful for dynamic vegetation models.

## Gustiness

In unstable conditions, especially when the mean wind speed approaches zero (free convection limit), convective eddies generated by surface heating maintain turbulent exchange. Standard bulk formulas using only the mean wind speed difference ($\Delta U$) would erroneously predict zero fluxes.

To account for this, `SurfaceFluxes.jl` uses an **effective wind speed** ($U_{\text{eff}}$):

```math
U_{\text{eff}} = \max\left( \left(\Delta u^2 + \Delta v^2\right)^{1/2}, \; U_{\text{gust}} \right)
```

where $U_{\text{gust}}$ is a parameterized gustiness velocity scale representing the contribution of sub-grid eddies.

### Parameterizations and Dispatch

The gustiness formulation is controlled by the `GustinessSpec` type in the `SurfaceFluxInputs`. This design allows for **type-stable dispatch** and straightforward broadcasting over heterogeneous surfaces. For example, a model coupled to both land and ocean can use an array of input structs where some elements use `ConstantGustinessSpec` (land) and others use `DeardorffGustinessSpec` (ocean). The solver simply calls `surface_fluxes` and correct method is dispatched automatically.

#### 1. Constant Gustiness

```julia
ConstantGustinessSpec(value)
```

Uses a fixed tuning parameter, e.g., $U_{\text{gust}} = 1.0 \, \mathrm{m~s^{-1}}$. Often used over land surfaces.

#### 2. Deardorff Gustiness

```julia
DeardorffGustinessSpec()
```

Scales $U_{\text{gust}}$ with the convective velocity scale $w_*$, following [Deardorff (1970)](https://doi.org/10.1175/1520-0469(1970)027<1211:CVATSF>2.0.CO;2) and [Beljaars (1995)](https://doi.org/10.1002/qj.49712152203). This is physically robust for the unstable boundary layer over fluid surfaces (ocean/lakes).

```math
U_{\text{gust}} = \beta w_* = \beta (B z_i)^{1/3}
```

where:

- The surface buoyancy flux is denoted by $B$.
- The boundary layer height is $z_i$.
- The scaling coefficient is $\beta$ (typically $\approx 1.0$).

Since $B$ depends on the fluxes, and the fluxes depend on $U_{\text{eff}}$ (and thus $B$), this introduces a nonlinear coupling that is resolved by an iterative solver.

## Reference

- Beljaars, A. C. M. (1995). The parametrization of surface fluxes in large-scale models under free convection. *Quarterly Journal of the Royal Meteorological Society*, 121, 255-270. [DOI: 10.1002/qj.49712152203](https://doi.org/10.1002/qj.49712152203)

- Deardorff, J. W. (1970). Convective velocity and temperature scales for the unstable planetary boundary layer and for Rayleigh convection. *Journal of the Atmospheric Sciences*, 27, 1211-1213. [DOI: 10.1175/1520-0469(1970)027<1211:CVATSF>2.0.CO;2](https://doi.org/10.1175/1520-0469(1970)027<1211:CVATSF>2.0.CO;2)

- Fairall, C. W., Bradley, E. F., Hare, J. E., Grachev, A. A., & Edson, J. B. (2003). Bulk parameterization of air–sea fluxes: Updates and verification for the COARE algorithm. *Journal of Climate*, 16, 571-591. [DOI: 10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2](https://doi.org/10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2)

- Raupach, M. R. (1994). Simplified expressions for vegetation roughness length and zero-plane displacement as functions of canopy height and area index. *Boundary-Layer Meteorology*, 71, 211-216. [DOI: 10.1007/BF00709229](https://doi.org/10.1007/BF00709229)

- Yatunin, D., Byrne, S., Kawczynski, C., Kandala, S., Bozzola, G., Sridhar, A., Shen, Z., Jaruga, A., Sloan, J., He, J., Huang, D.Z., Barra, V., Knoth, O., Ullrich, P., Schneider, T., 2026: The CliMA atmosphere dynamical core: Concepts, numerics, and scaling. *Journal of Advances in Modeling Earth Systems*, in press. [DOI:10.1029/2025MS005014](https://doi.org/10.1029/2025MS005014)
