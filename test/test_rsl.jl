using Test

import SurfaceFluxes
const SF = SurfaceFluxes
SurfaceFluxes.error_on_non_convergence() = true

import SurfaceFluxes.UniversalFunctions as UF
using StaticArrays
using Thermodynamics
Thermodynamics.print_warning() = false

include(joinpath(pkgdir(SurfaceFluxes), "parameters", "create_parameters.jl"))
const SFP = SF.Parameters
FT = Float64;

### Generate parameter lists
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
param_set = create_parameters(toml_dict, UF.Gryanik())
thermo_params = SFP.thermodynamics_params(param_set)
uft = SFP.universal_func_type(param_set)
### 

""" 
Parameters from Bonan (2019) Chapter 6

# References
- [Bonan2019](@cite) (Chapter 6)

# Original Research
- [Physick1995](@cite)

"""
T₀ = FT(273.15)
z_star = FT(49)
d = FT(19)
h_c = FT(22)
z0m = FT(0.6) # Figure 6.4
z0b = FT(0.135) * z0m
κ = SFP.von_karman_const(param_set)

u_sfc = FT(3.7)
u_in = FT(0)
θ_sfc = FT(289.7) # K
θ_in = FT(298.0) # K

ts_sfc_test = Thermodynamics.PhaseEquil_pθq(thermo_params, 99340.0, 289.7, 0.0107)
ts_int_test = Thermodynamics.PhaseEquil_pθq(thermo_params, 95342.0, 298.0, 0.0085)
state_in = SF.InteriorValues(FT(350), (FT(0.0), FT(0)), ts_int_test)
state_sfc = SF.SurfaceValues(FT(0), (FT(0), FT(0)), ts_sfc_test)
sc = SF.ValuesOnly{FT}(; state_in, state_sfc, z0m, z0b)
Z = collect(range(FT(d + 0.1), stop = FT(100), length = 250))


# TODO: expand test to include multiple canopies
@testset "De Ridder and Physick RSL Profiles Agree" begin
    testcanopy = SurfaceFluxes.SparseCanopy{FT}(d, z_star)
    PhysickRSL = SurfaceFluxes.PhysickRSL(testcanopy)
    DeRidderRSL = SurfaceFluxes.DeRidderRSL(testcanopy)

    @testset "Wind Profiles" begin
        for (iz, z) in enumerate(Z)
            for L_MO in FT[-1000, -50, -10, 30, 50, 1000]
                du_ridder = SF.recover_profile(
                    param_set,
                    sc,
                    L_MO,
                    FT(z),
                    u_in,
                    u_sfc,
                    UF.MomentumTransport(),
                    uft,
                    SF.FDScheme(),
                    DeRidderRSL,
                )
                du_physick = SF.recover_profile(
                    param_set,
                    sc,
                    L_MO,
                    FT(z),
                    u_in,
                    u_sfc,
                    UF.MomentumTransport(),
                    uft,
                    SF.FDScheme(),
                    PhysickRSL,
                )

                @test isapprox(du_ridder, du_physick, rtol = 0.25)
            end
        end
    end

    @testset "Temperature Profiles" begin
        for (iz, z) in enumerate(Z)
            for L_MO in FT[-1000, -50, -10, 30, 50, 1000]
                dθ_ridder = SF.recover_profile(
                    param_set,
                    sc,
                    L_MO,
                    FT(z),
                    θ_in,
                    θ_sfc,
                    UF.HeatTransport(),
                    uft,
                    SF.FDScheme(),
                    DeRidderRSL,
                )
                dθ_physick = SF.recover_profile(
                    param_set,
                    sc,
                    L_MO,
                    FT(z),
                    θ_in,
                    θ_sfc,
                    UF.HeatTransport(),
                    uft,
                    SF.FDScheme(),
                    PhysickRSL,
                )

                @test isapprox(dθ_ridder, dθ_physick, rtol = 0.25)
            end
        end
    end
end
