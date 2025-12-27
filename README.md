<div align="center">
  <img src="docs/src/assets/logo.svg" alt="SurfaceFluxes.jl Logo" width="128" height="128">
</div>

# SurfaceFluxes.jl

A package for computing surface fluxes between the atmosphere, ocean, and land using Monin-Obukhov Similarity Theory (MOST).

SurfaceFluxes.jl provides robust, efficient methods for calculating turbulent surface fluxes of momentum, heat, and moisture. It supports GPU broadcasting, automatic differentiation, and multiple universal function parameterizations (Businger, Gryanik, Grachev), making it ideal for high-performance climate modeling.

|                           |                                                                          |
|--------------------------:|:-------------------------------------------------------------------------|
| **Stable Release**        | [![stable][stable-img]][stable-url] [![docs-stable][docs-stable-img]][docs-stable-url] |
| **Latest Documentation**  | [![dev][docs-dev-img]][docs-dev-url]                                     |
| **Unit Tests**            | [![unit tests][gha-ci-img]][gha-ci-url] [![codecov][codecov-img]][codecov-url] |
| **Downloads**             | [![Downloads][dlt-img]][dlt-url]                                         |

[stable-img]: https://img.shields.io/github/v/release/CliMA/SurfaceFluxes.jl?label=stable
[stable-url]: https://github.com/CliMA/SurfaceFluxes.jl/releases/latest

[docs-stable-img]: https://img.shields.io/badge/docs-stable-green.svg
[docs-stable-url]: https://CliMA.github.io/SurfaceFluxes.jl/stable/

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://CliMA.github.io/SurfaceFluxes.jl/dev/

[gha-ci-img]: https://github.com/CliMA/SurfaceFluxes.jl/actions/workflows/ci.yml/badge.svg
[gha-ci-url]: https://github.com/CliMA/SurfaceFluxes.jl/actions/workflows/ci.yml

[codecov-img]: https://codecov.io/gh/CliMA/SurfaceFluxes.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/CliMA/SurfaceFluxes.jl

[dlt-img]: https://img.shields.io/badge/dynamic/json?url=http%3A%2F%2Fjuliapkgstats.com%2Fapi%2Fv1%2Ftotal_downloads%2FSurfaceFluxes&query=total_requests&label=Downloads
[dlt-url]: https://juliapkgstats.com/pkg/SurfaceFluxes

## Features

- **Monin-Obukhov Similarity Theory**: Robust iterative solver for stability-dependent surface fluxes
- **Universal Function Parameterizations**: [Businger et al. (1971)](https://doi.org/10.1175/1520-0469(1971)028%3C0181:FPRITA%3E2.0.CO;2), [Gryanik et al. (2020)](https://doi.org/10.1175/JAS-D-19-0255.1), and [Grachev et al. (2007)](https://doi.org/10.1007/s10546-007-9177-6) formulations
- **GPU Support**: Full GPU acceleration with CUDA.jl and other GPU array types
- **Land and Ocean Parameterizations**: Support for parameterizations for land and ocean surfaces, including roughness lengths that depend on wind speed (ocean) and vegetation characteristics (land)
- **Dynamic Skin States**: Supports dynamic calculations of skin temperatures and humidities via user-supplied functions
- **Finite-Difference and Finite-Volume Schemes**: Supports both finite-difference (point-wise) and finite-volume (layer-averaged) formulations following [Nishizawa & Kitamura (2018)](https://doi.org/10.1029/2018MS001534)
- **AD Compatible**: Works with automatic differentiation frameworks for integration into differentiable models

## Quick Example

```julia
using SurfaceFluxes
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions as UF
import ClimaParams 

# Create parameters
FT = Float64
param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)

# Compute surface fluxes
result = surface_fluxes(
    param_set,
    T_int,           # Interior temperature [K]
    q_tot,           # Interior total specific humidity [kg/kg]
    q_liq,           # Interior liquid specific humidity [kg/kg]
    q_ice,           # Interior ice specific humidity [kg/kg]
    ρ_int,           # Interior density [kg/m³]
    T_sfc,           # Surface temperature [K]
    q_sfc,           # Surface specific humidity [kg/kg]
    Φ_sfc,           # Surface geopotential [m²/s²]
    Δz,              # Height above surface [m]
    d,               # Displacement height [m]
    u_int,           # Interior wind (u, v) [m/s]
    u_sfc,           # Surface wind (u, v) [m/s]
)

# Access results
result.shf      # Sensible heat flux [W/m²]
result.lhf      # Latent heat flux [W/m²]
result.E        # Evaporation rate [kg/(m²·s)]
result.ustar    # Friction velocity [m/s]
result.ρτxz     # Momentum flux, x-component [N/m²]
result.ρτyz     # Momentum flux, y-component [N/m²]
result.Cd       # Drag coefficient
result.Ch       # Heat exchange coefficient
```
