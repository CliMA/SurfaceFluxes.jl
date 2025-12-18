# API Reference

## Main Solver Interface

```@docs
SurfaceFluxes.surface_fluxes
SurfaceFluxes.SurfaceFluxConditions
SurfaceFluxes.SurfaceFluxConfig
SurfaceFluxes.FluxSpecs
SurfaceFluxes.SurfaceFluxInputs
SurfaceFluxes.SolverOptions
SurfaceFluxes.SolverScheme
SurfaceFluxes.PointValueScheme
SurfaceFluxes.LayerAverageScheme
SurfaceFluxes.compute_profile_value
```

## Flux Calculations

Functions for computing specific fluxes.

```@docs
SurfaceFluxes.sensible_heat_flux
SurfaceFluxes.latent_heat_flux
SurfaceFluxes.buoyancy_flux
SurfaceFluxes.evaporation
SurfaceFluxes.momentum_fluxes
SurfaceFluxes.state_bulk_richardson_number
```

## Exchange Coefficients

Non-dimensional exchange coefficients and conductances.

```@docs
SurfaceFluxes.drag_coefficient
SurfaceFluxes.heat_exchange_coefficient
SurfaceFluxes.heat_conductance
```

## Physical Scales & Variances

Functions for computing Monin-Obukhov similarity scales and variances.

```@docs
SurfaceFluxes.compute_physical_scale_coeff
SurfaceFluxes.compute_ustar
SurfaceFluxes.compute_theta_star
SurfaceFluxes.compute_q_star
SurfaceFluxes.u_variance
SurfaceFluxes.scalar_variance
SurfaceFluxes.theta_variance
SurfaceFluxes.obukhov_length
SurfaceFluxes.obukhov_stability_parameter
```

## Roughness & Gustiness

```@docs
SurfaceFluxes.roughness_lengths
SurfaceFluxes.gustiness_constant
SurfaceFluxes.ConstantRoughnessParams
SurfaceFluxes.COARE3RoughnessParams
SurfaceFluxes.RaupachRoughnessParams
SurfaceFluxes.ConstantGustinessSpec
SurfaceFluxes.DeardorffGustinessSpec
SurfaceFluxes.MoistModel
SurfaceFluxes.DryModel
```

## Universal Functions

The `UniversalFunctions` sub-module defines the stability functions $\phi(\zeta)$ and $\psi(\zeta)$.

```@docs
SurfaceFluxes.UniversalFunctions
SurfaceFluxes.UniversalFunctions.phi
SurfaceFluxes.UniversalFunctions.psi
SurfaceFluxes.UniversalFunctions.Psi
```

### Parameter Types

```@docs
SurfaceFluxes.UniversalFunctions.BusingerParams
SurfaceFluxes.UniversalFunctions.GryanikParams
SurfaceFluxes.UniversalFunctions.GrachevParams
```

### Transport Types

```@docs
SurfaceFluxes.UniversalFunctions.MomentumTransport
SurfaceFluxes.UniversalFunctions.HeatTransport
```
