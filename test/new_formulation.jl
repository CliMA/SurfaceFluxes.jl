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

function z0test(surface_args, similarity_scales, atmos_state, param_set)
    u★ = similarity_scales.momentum
    FT =typeof(u★)
    return FT(0.015 * u★^2 / 9.81)
end


atmos_state = AtmosState(
                  FT(5),
                  FT(1),
                  FT(0.002),
                  FT(298),
                  FT(15),
                  FT(1), # gustiness needs to be a function of u,v, ustar
                  FT(100),
                  (ρ=FT(1.20), arg𝑏=FT(0.01), arg𝑐=z0test),
                )

surface_state = SurfaceState(
                  (𝑧0m=FT(0.01), 𝑧0θ=FT(0.01), 𝑧0q=z0test),
                  FT(0),
                  FT(0),
                  FT(0.003),
                  FT(304),
                  FT(0),
                  (arg𝑎=FT(0.01), arg𝑏=FT(0.01), arg𝑐=z0test),
                )


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

similarity_profile = ufunc
similarity_scales = refine_similarity_variables(Σ_est, ΔU_est, 
                                     similarity_profile,
                                     surface_state, 
                                     atmos_state, param_set)

fluxes = compute_similarity_theory_fluxes(similarity_profile, 
                                 surface_state,
                                 atmos_state, 
                                 param_set)

#### Diagnostics
@info ufunc
@info "With ufunc.L = $(ufunc.L) the Monin Obukhov length"
@info fluxes.sensible_heat
@info fluxes.latent_heat
@info fluxes.water_vapor
@info fluxes.x_momentum
@info fluxes.y_momentum
@info fluxes.r_ae
@info fluxes.scale_vars
