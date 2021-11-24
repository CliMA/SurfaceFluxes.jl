using NCDatasets
using Plots
using SurfaceFluxes
const SF = SurfaceFluxes
const UF = SurfaceFluxes.UniversalFunctions
using Statistics
using StaticArrays
using Thermodynamics
const TD = Thermodynamics

using CLIMAParameters: AbstractEarthParameterSet
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

files = ["DYCOMS_RF01.nc", "Bomex.nc", "Rico.nc", "Gabls.nc"];
for (n, f) in enumerate(files)
    @info "Casename = $f"
    data = NCDataset("./test/" * f)
    z = data["z_half"]
    u = data.group["profiles"]["u_mean"]
    v = data.group["profiles"]["v_mean"]
    qt = data.group["profiles"]["qt_mean"]
    ql = data.group["profiles"]["ql_mean"]
    ρ = data.group["profiles"]["rho"]
    p0 = data.group["profiles"]["p0"]
    b = data.group["profiles"]["buoyancy_mean"]
    θli = data.group["profiles"]["thetali_mean"] # Care with variable names across cases ()?
    T = data.group["profiles"]["temperature_mean"] # Care with variable names across cases ()?  
    SHF = data.group["timeseries"]["shf_surface_mean"] # Care with variable names across cases ()?
    LHF = data.group["timeseries"]["lhf_surface_mean"] # Care with variable names across cases ()?

    function getval(X)
        X[:, end]
    end

    # Data at first interior node (x_in)
    # TODO Make sure that the first node ii = 1 is at the "surface" 
    # for the tests to be consistent
    ii = 2
    FT = Float32
    ArrayType = Array

    shf = mean(getval(SHF))
    lhf = mean(getval(LHF))

    z_in = Tuple(z)[ii]
    u_in = getval(u)[ii]
    v_in = getval(v)[ii]
    b_in = Tuple(getval(b))[ii]
    ρ_in = Tuple(getval(ρ))[ii]
    θ_in = Tuple(getval(θli))[ii]
    T_in = Tuple(getval(T))[ii]
    qt_in = Tuple(getval(qt))[ii]

    ## Initial guesses for MO parameters
    LMO_init = eps(FT)
    u_star_init = FT(0.1)
    b_star_init = -FT(0.1)
    qt_star_init = -FT(1e-5)

    # Surface values for variables
    jj = 1
    u_sfc = FT(0)
    v_sfc = FT(0)
    b_sfc = Tuple(getval(b))[jj]
    ρ_sfc = Tuple(getval(ρ))[jj]
    qt_sfc = Tuple(getval(qt))[jj]
    θ_sfc = Tuple(getval(θli))[jj]
    T_sfc = Tuple(getval(T))[jj]
    z_sfc = FT(0)

    # Roughness
    ## Roughness lengths
    z0m = FT(0.001)
    z0b = FT(0.001)

    ts_sfc = TD.PhaseEquil_ρθq(param_set, ρ_sfc, θ_sfc, qt_sfc)
    ts_in = TD.PhaseEquil_ρθq(param_set, ρ_in, θ_in, qt_in)

    u_in = SVector{2, FT}(u_in, v_in)
    u_sfc = SVector{2, FT}(0, 0)

    state_gabls_sfc = SF.ValuesSurface{FT, SVector, TD.ThermodynamicState}(
        ts_sfc,
        u_sfc,
        z_sfc,
    )
    state_gabls_in =
        SF.ValuesInterior{FT, SVector, TD.ThermodynamicState}(ts_in, u_in, z_in)

    state_rico_sfc = SF.ValuesSurface{FT, SVector, TD.ThermodynamicState}(
        ts_sfc,
        u_sfc,
        z_sfc,
    )
    state_rico_in =
        SF.ValuesInterior{FT, SVector, TD.ThermodynamicState}(ts_in, u_in, z_in)

    state_bomex_sfc = SF.ValuesSurface{FT, SVector, TD.ThermodynamicState}(
        ts_sfc,
        u_sfc,
        z_sfc,
    )
    state_bomex_in =
        SF.ValuesInterior{FT, SVector, TD.ThermodynamicState}(ts_in, u_in, z_in)

    state_dycoms_sfc = SF.ValuesSurface{FT, SVector, TD.ThermodynamicState}(
        ts_sfc,
        u_sfc,
        z_sfc,
    )
    state_dycoms_in =
        SF.ValuesInterior{FT, SVector, TD.ThermodynamicState}(ts_in, u_in, z_in)

    # shf and lhf in dycoms are mising /wrong in the file
    if f == "DYCOMS_RF01.nc"
        sc = SF.GivenFluxes(
            state_dycoms_in,
            state_dycoms_sfc,
            FT(15),
            FT(115),
            UF.Businger,
            z0m,
            z0b,
        )
    elseif f == "Bomex.nc"
        sc = SF.GivenFluxesAndFrictionVelocity(
            state_bomex_in,
            state_bomex_sfc,
            shf,
            lhf,
            FT(0.28),
            UF.Businger,
            z0m,
            z0b,
        )
    elseif f == "Rico.nc"
        sc = SF.GivenCoefficients(
            state_rico_in,
            state_rico_sfc,
            FT(0.001229),
            FT(0.001094),
            UF.Businger,
            z0m,
            z0b,
        )
    elseif f == "Gabls.nc"
        sc = SF.GivenSurfaceValues(
            state_gabls_in,
            state_gabls_sfc,
            UF.Businger,
            z0m,
            z0b,
            LMO_init,
        )
    end
    result = surface_conditions(param_set, sc)
end
