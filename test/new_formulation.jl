include("../src/similarity_theory.jl")

using Test
import Thermodynamics as TD
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions.BusingerParams
import SurfaceFluxes.UniversalFunctions as UF
import ClimaParams as CP

# Assume constant values
# For method tests define functions for roughness that return constant values


FT = Float32 # TODO Check all floattypes
Σ₀ = SimilarityScales{FT, FT, FT}(1e-4,1e-4,1e-4)
Σₜ = SimilarityScales{FT, FT, FT}(1e-5,1e-5,1e-5)
ΔΣ = Σₜ- Σ₀

# Parameters (ClimaParams types)
param_set = SFP.SurfaceFluxesParameters(FT, BusingerParams)
thermo_params = param_set.thermo_params

# Assign states

z0test(u★, ζ) = FT(0.01)   

atmos_state = AtmosState(
                  FT(1),
                  FT(1),
                  FT(0.002),
                  FT(298),
                  FT(15),
                  FT(1), # gustiness needs to be a function of u,v, ustar
                  FT(100),
                  (arg𝑎=FT(0.01), arg𝑏=FT(0.01), arg𝑐=z0test),
                )

surface_state = SurfaceState(
                  (𝑧0m=FT(0.01), 𝑧0θ=FT(0.01), 𝑧0q=z0test),
                  FT(0),
                  FT(0),
                  FT(0.002),
                  FT(299),
                  FT(0),
                  (arg𝑎=FT(0.01), arg𝑏=FT(0.01), arg𝑐=z0test),
                )


# Test function inputs within `args` : see ClimaOcean for uniformity 
# in unpack methods.
@assert atmos_state.args.arg𝑐(1,2) == FT(0.01)


# Line by line debug and test for `refine` function
gustiness = atmos_state.gustiness_parameter
Δstate = state_differences(surface_state, atmos_state, Σ₀, param_set); 
(; 𝑧0m, 𝑧0θ, 𝑧0q) = surface_state.roughness_lengths

# Generic info block 

ζ₀ = FT(-10)
L★ = Δstate.Δh ./ ζ₀
sfc_params = SFP.uf_params(param_set)
similarity_theory = SFP.universal_func_type(param_set)
ufunc = UF.universal_func(similarity_theory, L★, sfc_params)
# We shouldn't need both these! ζ and `BusingerParams` should be enough (at the user level)
# to define all similarity functions 𝜓, 𝜙

Σ_est = (momentum=FT(0.1),temperature=FT(0.01),water_vapor=FT(0.001))
ΔU_est = FT(10)

# Initial guess
u★ = Σ_est.momentum
θ★ = Σ_est.temperature
q★ = Σ_est.water_vapor
uτ = ΔU_est
## TODO Define methods for buoyancy_scale; current implementation uses `compute_bstar`
#b★ = buoyancy_scale(θ★, q★, thermo_params)
b★ = FT(0.2)

## TODO Fix Parameter unpack methods (unify between ClimaOcean and ClimaParams)
𝑔 = FT(9.81)
𝜅 = FT(0.4)
L★ = ifelse(b★ == 0, zero(b★), - u★^3 * atmos_state.θ_a / (u★ * θ★ * 𝜅 * 𝑔))
ζ = Δstate.Δh / L★ 
ψm = UF.psi(ufunc, ζ, UF.MomentumTransport())
ψs = UF.psi(ufunc, ζ, UF.HeatTransport()) # TODO Rename HeatTransport > ScalarTransport
ψm₀ = UF.psi(ufunc, 𝑧0m * ζ / Δstate.Δh, UF.MomentumTransport())
ψh₀ = UF.psi(ufunc, 𝑧0θ * ζ / Δstate.Δh, UF.HeatTransport())
ψq₀ = UF.psi(ufunc, 𝑧0q(u★,ζ) * ζ / Δstate.Δh, UF.HeatTransport())

# compute rhs in Δχ/u★ = (f(ζ,𝑧0...))
F_m = log(Δstate.Δh / 𝑧0m) - ψm + ψm₀
F_h = log(Δstate.Δh / 𝑧0θ) - ψs + ψh₀
F_q = log(Δstate.Δh / 𝑧0q(u★, ζ)) - ψs + ψq₀

# Review against nishizawa notation
χu = 𝜅/F_m 
χθ = 𝜅/F_m
χq = 𝜅/F_m

u★ = χu * uτ
θ★ = χθ * Δstate.Δθ
q★ = χq * Δstate.Δq

# Buoyancy flux similarity scale for gustiness (Edson 2013)
h_atmos_boundary_layer = FT(100)
hᵢ = h_atmos_boundary_layer
Jᵇ = - u★ * b★
Uᴳ = gustiness * cbrt(Jᵇ * hᵢ)

# New velocity difference accounting for gustiness
ΔU = sqrt(Δstate.Δu^2 + Δstate.Δv^2 + Uᴳ^2)

# TODO: z0test to be redefined with `surface_args`, `similarity_scales` as args

similarity_profile = ufunc
compute_similarity_theory_fluxes(similarity_profile, 
                                 surface_state,
                                 atmos_state, 
                                 param_set)

#### Diagnostics
@info atmos_state.args
@info Δstate
@info propertynames(surface_state)
@info propertynames(atmos_state)
@info similarity_theory
@info ufunc
@info "With ufunc.L = $(ufunc.L) the Monin Obukhov length"
