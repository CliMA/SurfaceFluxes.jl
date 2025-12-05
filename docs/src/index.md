# SurfaceFluxes.jl

See the [API Reference](@ref) for complete API documentation.

## Source Code Structure

The `src/` directory contains the following files organized by functionality:

### Core Module Files
- **`SurfaceFluxes.jl`**: Main module file containing the primitive `surface_fluxes` API and the `obukhov_similarity_solution` solver
- **`types.jl`**: Type definitions including solver options, `SurfaceFluxConfig`, and iteration-state records
- **`utilities.jl`**: Helper functions for Richardson number computation and thermodynamic differences
- **`Parameters.jl`**: Parameter set definitions for physical constants
- **`UniversalFunctions.jl`**: Universal function implementations (Businger, Gryanik, Grachev)

### Flux and Exchange Coefficient Methods
- **`physical_scale_coefficient_methods.jl`**: Computation of physical scale coefficients for finite difference (Byun 1990) and finite volume (Nishizawa 2018) schemes
- **`friction_velocity_methods.jl`**: Friction velocity (u★) computation methods
- **`exchange_coefficients.jl`**: Heat conductance, momentum exchange coefficient (Cd), and heat exchange coefficient (Ch) computations
- **`bulk_fluxes.jl`**: Sensible heat flux, latent heat flux, evaporation rate, buoyancy flux, and momentum flux computations

### Surface Configuration Helpers
- **`input_builders.jl`**: Normalization of primitive inputs plus parameterization selection via `SurfaceFluxConfig`
- **`roughness_lengths.jl`**: Roughness length evaluation for built-in parameterizations

### Profile Recovery
- **`profile_recovery.jl`**: Functions to recover vertical profiles within the surface layer using Monin-Obukhov similarity theory

For detailed API documentation of all types and methods, see the [API Reference](@ref).

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

## Documentation
- [API Reference](API.md)
- [Universal Functions](UniversalFunctions.md)
