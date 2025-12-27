# Test Suite

*Comprehensive testing for correctness and numerical stability*

## Overview

SurfaceFluxes.jl employs a multi-tiered testing strategy to ensure physical correctness, numerical robustness, and software compatibility (GPU/AD). The suite covers:
1.  **Regression Tests:** Preventing breaking changes against known baselines.
2.  **Universal Functions Tests:** Verifying mathematical properties of stability corrections.
3.  **Physical Consistency:** Checking flux directions and energy balance.
4.  **Software Compatibility:** Ensuring compatibility with GPU execution (CUDA.jl) and Automatic Differentiation (ForwardDiff.jl).

## Automatic Differentiation (AD)

The package is designed to be fully differentiable. This is critical for coupled models where surface fluxes are part of a larger implicit solver or optimization loop.

### AD Test Example

We verify AD compatibility by comparing the derivatives computed via `ForwardDiff.jl` against finite-difference approximations.

```julia
using ForwardDiff
using SurfaceFluxes

# Wrapper function for the flux calculation
function compute_shf(T_sfc)
    # create inputs with T_sfc as a Dual number...
    result = surface_fluxes(param_set, ..., T_sfc, ...)
    return result.shf
end

# Compute derivative d(SHF)/d(T_sfc)
dSHF_dT = ForwardDiff.derivative(compute_shf, 300.0)
```

Tests ensure that:
1.  Dual numbers propagate through the entire solver (no type conversion errors).
2.  Derivatives are accurate compared to finite differences across stable and unstable regimes.

## Detailed Test Categories

### Regression Tests

Predefined test cases with known expected outputs, ensuring that code changes do not introduce regressions. Covers over 1600 cases across different stability regimes, floating-point types, and roughness parameterizations.

These tests are also valuable for performance tuning of hyperparameters, such as determining the minimum `maxiter` required to achieve a target accuracy (e.g., < 10% error in fluxes).

### Universal Functions Tests

Rigorous verification of the stability correction functions:

#### Type Stability
**Test:** `Type stability`

Verifies that all functions return values of the correct floating-point type (`Float32` or `Float64`) matching the input parameter type.

#### Neutral Limit Behavior
**Test:** `Neutral logarithmic velocity profile`

Verifies that in the neutral limit ($L \to \infty$, i.e., $\zeta \to 0$), the velocity profile collapses to the logarithmic law of the wall:
```math
\begin{equation}
u(z) = \frac{u_*}{\kappa} \ln\left(\frac{z-d}{z_0}\right)
\end{equation}
```

#### Asymptotic Behavior
**Test:** `Asymptotic behavior (|ζ| → ∞)`

For very stable conditions ($\zeta \gg 1$), the functions should approach their asymptotic limits:
- **Gryanik**: $\phi_m(\zeta) \sim \zeta^{1/3}$
- **Grachev**: $\phi_m(\zeta) \sim \zeta^{1/3}$

#### Mathematical Consistency
Several tests ensure the internal consistency of the definitions:
- **Derivative Consistency**: Checks $\phi(\zeta) \approx \phi(0) - \zeta \cdot \psi'(\zeta)$ using finite differences.
- **Integral Consistency**: Checks $\psi(\zeta) \approx \int (\phi(0) - \phi(x))/x \, dx$ using numerical quadrature.
- **Continuity**: Verifies continuous transitions in functions and their derivatives, **especially at the neutral limit** ($\zeta \to 0$).

## Running the Tests

```julia
using Pkg
Pkg.test("SurfaceFluxes")
```

Or run specific test files:

```julia
julia --project=test test/test_universal_functions.jl
```

## Quality Assurance

The package uses [Aqua.jl](https://github.com/JuliaTesting/Aqua.jl) for automated quality checks including ambiguity detection and unbound type parameters.
