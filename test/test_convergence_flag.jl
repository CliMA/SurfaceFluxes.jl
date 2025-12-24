using Test
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD
import Thermodynamics.Parameters as TP
import ClimaParams as CP

FT = Float32

# Create parameters
toml_dict = CP.create_toml_dict(FT)
param_set = SFP.SurfaceFluxesParameters(toml_dict, UF.BusingerParams)
thermo_params = SFP.thermodynamics_params(param_set)

@testset "Convergence Flag Verification" begin
    # 1. Standard Case - Should Converge
    T_int = FT(300)
    T_sfc = FT(301)
    q_tot_int = FT(0.01)
    q_vap_sfc = FT(0.012)
    ρ_int = FT(1.0)
    u_int = (FT(5), FT(0))
    u_sfc = (FT(0), FT(0))
    Φ_sfc = FT(0)
    Δz = FT(10)
    d = FT(0)

    # Standard solver options
    opts = SF.SolverOptions{FT}(maxiter = 10, tol = 1e-4, forced_fixed_iters = false)


    # Defaults for positional args
    roughness = nothing
    config = SF.default_surface_flux_config(FT)
    scheme = SF.PointValueScheme()

    sf = SF.surface_fluxes(
        param_set,
        T_int, q_tot_int, FT(0), FT(0), ρ_int,
        T_sfc, q_vap_sfc,
        Φ_sfc, Δz, d,
        u_int, u_sfc,
        roughness,
        config,
        scheme,
        opts,
    )

    @test sf.converged == true
    @test sf.L_MO != 0
    @test sf.shf != 0

    # 2. Forced Non-Convergence Case
    # By setting maxiter to 1 and starting far from the root, we expect it to fail convergence check
    opts_fail = SF.SolverOptions{FT}(maxiter = 1, tol = 1e-10, forced_fixed_iters = false)

    sf_fail = SF.surface_fluxes(
        param_set,
        T_int, q_tot_int, FT(0), FT(0), ρ_int,
        T_sfc, q_vap_sfc,
        Φ_sfc, Δz, d,
        u_int, u_sfc,
        roughness,
        config,
        scheme,
        opts_fail,
    )

    @test sf_fail.converged == false

    # 3. Prescribed Fluxes - Should be true
    flux_specs_prescribed = SF.FluxSpecs(ustar = FT(0.5), shf = FT(100.0), lhf = FT(200.0))
    sf_prescribed = SF.surface_fluxes(
        param_set,
        T_int, q_tot_int, FT(0), FT(0), ρ_int,
        T_sfc, q_vap_sfc,
        Φ_sfc, Δz, d,
        u_int, u_sfc,
        roughness,
        config,
        scheme,
        nothing, # solver_opts
        flux_specs_prescribed,
    )
    @test sf_prescribed.converged == true

    # 4. Prescribed Coefficients - Should be true
    flux_specs_coeffs = SF.FluxSpecs(Cd = FT(0.001), Ch = FT(0.001))
    sf_coeffs = SF.surface_fluxes(
        param_set,
        T_int, q_tot_int, FT(0), FT(0), ρ_int,
        T_sfc, q_vap_sfc,
        Φ_sfc, Δz, d,
        u_int, u_sfc,
        roughness,
        config,
        scheme,
        nothing, # solver_opts
        flux_specs_coeffs,
    )
    @test sf_coeffs.converged == true

end
