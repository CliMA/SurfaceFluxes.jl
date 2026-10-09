# API Reference

## Main Solver Interface

```@docs
SurfaceFluxes.surface_fluxes
SurfaceFluxes.SurfaceFluxConditions
SurfaceFluxes.SurfaceFluxConfig
SurfaceFluxes.FluxSpecs
SurfaceFluxes.SolverOptions
SurfaceFluxes.SolverScheme
SurfaceFluxes.PointValueScheme
SurfaceFluxes.LayerAverageScheme
SurfaceFluxes.compute_profile_value
SurfaceFluxes.screen_level_values
SurfaceFluxes.dimensionless_profile_value
```

## Inputs Container

Many internal functions operate on a normalized "inputs container" (a `NamedTuple`)
built from the user-facing arguments by
[`build_surface_flux_inputs`](@ref SurfaceFluxes.build_surface_flux_inputs).

```@docs
SurfaceFluxes.build_surface_flux_inputs
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
SurfaceFluxes.surface_tke
SurfaceFluxes.scalar_variance
SurfaceFluxes.theta_variance
SurfaceFluxes.obukhov_length
SurfaceFluxes.obukhov_stability_parameter
```

## Utilities

```@docs
SurfaceFluxes.surface_density
SurfaceFluxes.surface_geopotential
SurfaceFluxes.interior_geopotential
SurfaceFluxes.effective_height
SurfaceFluxes.ReferenceAboveSurface
SurfaceFluxes.ReferenceAboveApparentSink
SurfaceFluxes.reference_above_surface
SurfaceFluxes.interior_vapor_specific_humidity
SurfaceFluxes.reference_height_valid
SurfaceFluxes.check_reference_height
SurfaceFluxes.invalidate_unless
```

## Roughness & Gustiness

```@docs
SurfaceFluxes.ConstantRoughnessParams
SurfaceFluxes.COARE3RoughnessParams
SurfaceFluxes.RaupachRoughnessParams
SurfaceFluxes.momentum_roughness
SurfaceFluxes.displacement_height
SurfaceFluxes.frontal_area_index
SurfaceFluxes.canopy_area_index
SurfaceFluxes.raupach_displacement_fraction
SurfaceFluxes.raupach_roughness_fraction
SurfaceFluxes.ConstantGustinessSpec
SurfaceFluxes.DeardorffGustinessSpec
SurfaceFluxes.FlooredDeardorffGustinessSpec
SurfaceFluxes.MoistModel
SurfaceFluxes.DryModel
SurfaceFluxes.gustiness_value
SurfaceFluxes.minimum_wind_speed
SurfaceFluxes.without_floor
SurfaceFluxes.free_convection_wind_speed
SurfaceFluxes.virtual_pottemps
SurfaceFluxes.depends_on_ustar
SurfaceFluxes.compute_ustar_and_roughness
```

## Roughness Sublayer

Models for the roughness sublayer (RSL) correction, which accounts for the enhanced
turbulent mixing above tall roughness elements (plant and urban canopies). Pass the chosen
model via `SurfaceFluxConfig(roughness, gustiness, moisture_model, rsl_model)`.

```@docs
SurfaceFluxes.NoRoughnessSubLayer
SurfaceFluxes.LinearRSL
SurfaceFluxes.ExponentialRSL
SurfaceFluxes.rsl_corrected_profile
SurfaceFluxes.rsl_profile_correction
```

## Stability Cap

Caps on the stability parameter in stable conditions, which hold the exchange
coefficients at their values at the cap for more stable conditions. Pass the chosen
cap via `SurfaceFluxConfig(roughness, gustiness, moisture_model, rsl_model, stability_cap)`.

```@docs
SurfaceFluxes.NoStabilityCap
SurfaceFluxes.ConstantStabilityCap
SurfaceFluxes.MaxHeatFluxStabilityCap
SurfaceFluxes.max_heat_flux_stability
SurfaceFluxes.neutral_momentum_roughness
SurfaceFluxes.stability_cap_value
SurfaceFluxes.with_stability_cap
SurfaceFluxes.resolved_stability_cap
SurfaceFluxes.capped_stability
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
