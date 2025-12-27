# SurfaceFluxes.jl

*A package for computing surface fluxes using Monin-Obukhov Similarity Theory (MOST)*

## Overview

SurfaceFluxes.jl provides a robust and efficient set of tools for calculating turbulent surface fluxes of momentum, heat, and moisture between the surface and the atmosphere. It is designed for use in climate models, large-eddy simulations (LES), and offline analysis.

The package implements **Monin-Obukhov Similarity Theory (MOST)** to relate surface fluxes to gradients of mean variables (wind speed, temperature, humidity) in the surface layer. It solves the nonlinear equations iteratively to determine the Obukhov stability parameter $\zeta$ and the corresponding fluxes.

### Key Features

- **Robust Iterative Solver**: Efficiently solves for the Monin-Obukhov stability parameter $\zeta$, supporting various roughness length parameterizations and flexible computation of skin temperature and humidity. 
- **Universal Functions**: Supports multiple parameterizations:
    - **Businger**: The classic Businger-Dyer formulations ([Businger et al. 1971](https://doi.org/10.1175/1520-0469(1971)028<0181:FPRITA>2.0.CO;2), [Dyer 1974](https://doi.org/10.1007/BF00240838)).
    - **Gryanik**: Improved functions for the stable boundary layer ([Gryanik et al. 2020](https://doi.org/10.1175/JAS-D-19-0255.1)).
    - **Grachev**: Functions derived from the SHEBA experiment for stable conditions over sea ice ([Grachev et al. 2007](https://doi.org/10.1007/s10546-007-9177-6)).
- **Thermodynamic Consistency**: Integrated with [Thermodynamics.jl](https://github.com/CliMA/Thermodynamics.jl) for accurate and consistent handling of moist air properties.
- **GPU Compatibility**: Type stable and designed for high-performance computing with full GPU support via [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl).
- **Automatic Differentiation**: Compatible with AD frameworks such as [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl).
- **Flexible Discretization**: Supports both **finite-difference** (point-wise) and **finite-volume** (layer-averaged) schemes.

## Installation

SurfaceFluxes.jl is a registered Julia package. You can install it using the package manager:

```julia
using Pkg
Pkg.add("SurfaceFluxes")
```

## Quick Start

Here is a simple example of how to compute surface fluxes given atmospheric state variables.

```@example quickstart
using SurfaceFluxes
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions as UF
import ClimaParams
import Thermodynamics as TD

# 1. Define Parameters
FT = Float32
# Create parameter set using Businger universal functions
param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)

# 2. Define State Variables
# Interior (air) state at height Δz
T_int = FT(298.0)      # Temperature [K]
q_tot = FT(0.017)      # Total specific humidity [kg/kg]
q_liq = FT(0.0)        # Liquid specific humidity [kg/kg]
q_ice = FT(0.0)        # Ice specific humidity [kg/kg]
ρ_int = FT(1.2)        # Air density [kg/m³]
u_int = (FT(5.0), FT(0.0)) # Wind vector [m/s]
Δz    = FT(25.0)       # Height above surface [m]

# Surface state
T_sfc = FT(300.0)      # Surface temperature [K]
q_sfc = FT(0.02)      # Surface specific humidity [kg/kg]
Φ_sfc = FT(0.0)        # Surface geopotential [m²/s²]
u_sfc = (FT(0.0), FT(0.0)) # Surface wind (usually 0)
d     = FT(5.0)        # Displacement height [m]

# 3. Compute Surface Fluxes
# This function iterates to find the stability parameter ζ
result = surface_fluxes(
    param_set,
    T_int, q_tot, q_liq, q_ice, ρ_int,
    T_sfc, q_sfc, Φ_sfc,
    Δz, d,
    u_int, u_sfc
)

# 4. Access Results
println("Sensible Heat Flux: ", round(result.shf; digits=1), " W/m²")
println("Latent Heat Flux:   ", round(result.lhf; digits=1), " W/m²")
println("Friction Velocity:  ", round(result.ustar; digits=3), " m/s")
println("Stability (ζ):      ", round(result.ζ; digits=4))
println("Obukhov Length:     ", round(result.L_MO; digits=1), " m")
```

The output shows positive sensible and latent heat fluxes (surface warmer and moister than air above it → upward sensible and latent heat transport). The stability parameter $\zeta < 0$ indicates unstable conditions.

## Automatic Differentiation

SurfaceFluxes.jl is fully differentiable, making it suitable for use in optimization loops or implicit solvers where Jacobians are required. It is tested with [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl).

Here is an example of computing the derivative of the sensible heat flux (SHF) with respect to surface temperature ($T_{\text{sfc}}$).

```@example quickstart
using ForwardDiff

# Define a wrapper function that takes T_sfc as the variable
function compute_shf(T_sfc_val)
    # Re-use parameters and other state variables from the Quick Start above
    # Note: ForwardDiff passes a dual number, so T_sfc_val will be of type Dual
    result = surface_fluxes(
        param_set,
        T_int, q_tot, q_liq, q_ice, ρ_int,
        T_sfc_val, q_sfc, Φ_sfc,
        Δz, d,
        u_int, u_sfc
    )
    return result.shf
end

# Compute the value and derivative at T_sfc = 300.0 K
T_sfc_val = FT(300.0)
shf_val = compute_shf(T_sfc_val)
dSHF_dT = ForwardDiff.derivative(compute_shf, T_sfc_val)

println("SHF at 300K:        ", round(shf_val; digits=2), " W/m²")
println("d(SHF)/d(Ts):       ", round(dSHF_dT; digits=2), " W/m²/K")
```

The derivative $\partial \text{SHF}/\partial T_s > 0$ indicates that increasing the surface temperature increases the upward sensible heat flux, as expected.

## Documentation Structure

| Section | Description |
|---------|-------------|
| [Surface Fluxes Theory](SurfaceFluxes.md) | Mathematical formulation of MOST |
| [Universal Functions](UniversalFunctions.md) | Stability functions (ϕ, ψ, Ψ) and their parameterizations |
| [Physical Scales](PhysicalScales.md) | Friction velocity $u_*$, temperature scale $\theta_*$, Obukhov length |
| [Exchange Fluxes](ExchangeFluxes.md) | Sensible heat, latent heat, and momentum fluxes |
| [Prescribed Conditions](PrescribedConditions.md) | Operating modes for prescribed fluxes, temperatures, or coefficients |
| [Test Suite](TestSuite.md) | Comprehensive testing for correctness and numerical stability |
| [API Reference](API.md) | Complete function documentation |

## Related Packages

SurfaceFluxes.jl is part of the [CliMA](https://github.com/CliMA) ecosystem:

- [Thermodynamics.jl](https://github.com/CliMA/Thermodynamics.jl): Moist thermodynamic calculations
- [RootSolvers.jl](https://github.com/CliMA/RootSolvers.jl): Iterative root-finding algorithms
- [ClimaParams](https://github.com/CliMA/ClimaParams): Centralized parameter management
