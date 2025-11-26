# SurfaceFluxes.jl

```@docs
SurfaceFluxes
```

## Source Code Structure

The `src/` directory contains the following files organized by functionality:

### Core Module Files
- **`SurfaceFluxes.jl`**: Main module file containing the `surface_fluxes` function and `obukhov_similarity_solution` solver
- **`types.jl`**: Type definitions including `SurfaceFluxInputs`, `SolverScheme`, and quantity/solver helper structs
- **`utilities.jl`**: Helper functions for state accessors, Richardson number computation, and thermodynamic differences
- **`Parameters.jl`**: Parameter set definitions for physical constants
- **`UniversalFunctions.jl`**: Universal function implementations (Businger, Gryanik, Grachev)

### Flux and Exchange Coefficient Methods
- **`physical_scale_coefficient_methods.jl`**: Computation of physical scale coefficients for finite difference (Byun 1990) and finite volume (Nishizawa 2018) schemes
- **`momentum_exchange_coefficient_methods.jl`**: Computation of momentum exchange coefficient (Cd) for neutral and stratified conditions
- **`heat_exchange_coefficient_methods.jl`**: Computation of heat exchange coefficient (Ch) for neutral and stratified conditions
- **`friction_velocity_methods.jl`**: Friction velocity (u★) computation methods
- **`sensible_heat_methods.jl`**: Sensible heat flux computations
- **`latent_heat_methods.jl`**: Latent heat flux computations
- **`buoyancy_flux_methods.jl`**: Buoyancy flux computations
- **`evaporation_methods.jl`**: Evaporation rate computations

### Surface Condition Input Types
- **`coefficient_inputs.jl`**: Methods for surface conditions specified via exchange coefficients (Cd, Ch)
- **`roughness_lengths.jl`**: Roughness length evaluation helpers for momentum/scalar callables

### Profile Recovery
- **`profile_recovery.jl`**: Functions to recover vertical profiles within the surface layer using Monin-Obukhov similarity theory

## Core input types

```@docs
SurfaceFluxes.StateValues
```

## Dispatch types

```@docs
SurfaceFluxes.Fluxes
SurfaceFluxes.FluxesAndFrictionVelocity
SurfaceFluxes.Coefficients
SurfaceFluxes.ValuesOnly
```

## User-facing methods

```@docs
SurfaceFluxes.surface_fluxes
SurfaceFluxes.recover_profile
```

# Parameters
Convenience constructors are provided for the `SurfaceFluxesParameters` and the various `UniversalFunctions` parameter structs.
To use them, you must first import ClimaParams:
```julia
import ClimaParams as CP
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions as UF

FT = Float64

# SurfaceFluxesParameters requires a float type and a UniversalFunctionsParameters type
SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)

# Or a TOML dict instead of a float type
toml_dict = CP.create_toml_dict(Float64)
SFP.SurfaceFluxesParameters(toml_dict, UF.GrachevParams)

# UniversalFunctionsParameters only require a float type or a TOML dict.
UF.BusingerParams(FT)
UF.GryanikParams(FT)
UF.GrachevParams(FT)
```

## Universal Functions

```@docs
SurfaceFluxes.UniversalFunctions
```

```@docs
SurfaceFluxes.UniversalFunctions.GryanikParams
SurfaceFluxes.UniversalFunctions.GrachevParams
SurfaceFluxes.UniversalFunctions.BusingerParams
```

```@docs
SurfaceFluxes.UniversalFunctions.phi
SurfaceFluxes.UniversalFunctions.psi
SurfaceFluxes.UniversalFunctions.Psi
```
