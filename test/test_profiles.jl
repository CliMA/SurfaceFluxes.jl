using NCDatasets
using Plots
using LaTeXStrings
using SurfaceFluxes
const SF = SurfaceFluxes
using Statistics
using Thermodynamics
const TD = Thermodynamics

using CLIMAParameters: AbstractEarthParameterSet
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

files = ["DYCOMS_RF01.nc", "Bomex.nc", "Rico.nc", "Gabls.nc"];
for (n, f) in enumerate(files)

    @info "Casename = $f"
    data = NCDataset("./test/"*f);
    z = data["z_half"];
    u = data.group["profiles"]["u_mean"];
    v = data.group["profiles"]["v_mean"];
    qt = data.group["profiles"]["qt_mean"];
    ql = data.group["profiles"]["ql_mean"];
    ρ = data.group["profiles"]["rho"];
    p0 = data.group["profiles"]["p0"];
    b = data.group["profiles"]["buoyancy_mean"];
    θli = data.group["profiles"]["thetali_mean"]; # Care with variable names across cases ()?
    T = data.group["profiles"]["temperature_mean"]; # Care with variable names across cases ()?

    function ave(X)
        return mean(X, dims=2)
    end

    # Data at first interior node (x_ave)
    # TODO Make sure that the first node ii = 1 is at the "surface" 
    # for the tests to be consistent
    ii = 2
    FT = Float32
    ArrayType = Array

    z_ave = Tuple(z)[ii]
    u_ave = Tuple(ave(u))[ii] 
    b_ave = Tuple(ave(b))[ii]
    ρ_ave = Tuple(ave(ρ))[ii]
    θ_ave = Tuple(ave(θli))[ii]
    T_ave = Tuple(ave(T))[ii]
    u_ave = Tuple(ave(u))[ii]
    qt_ave = Tuple(ave(qt))[ii]
    x_ave = ArrayType(FT[u_ave, b_ave])

    ## Initial guesses for MO parameters
    LMO_init = eps(FT)
    u_star_init = FT(0.1)
    b_star_init = -FT(0.1)
    qt_star_init = -FT(1e-5)
    MO_param_guess =
        ArrayType(FT[LMO_init, u_star_init, b_star_init, qt_star_init])

    # Surface values for variables
    jj = 1
    u_sfc = FT(0)
    b_sfc = Tuple(ave(b))[jj]
    ρ_sfc = Tuple(ave(ρ))[jj]
    qt_sfc = Tuple(ave(qt))[jj]
    θ_sfc = Tuple(ave(θli))[jj]
    T_sfc = Tuple(ave(T))[jj]
    z_sfc = FT(0)
    x_s = ArrayType(FT[u_sfc, b_sfc, qt_sfc])

    # Roughness
    ## Roughness lengths
    z0 = ArrayType(FT[
        5.86144925739178e-05,
        0.0001,
        0.000641655193293549,
        3.23383768877187e-05,
    ])
    zt = ArrayType(FT[
        3.69403636275411e-05,
        0.0001,
        1.01735489109205e-05,
        7.63933834969505e-05,
    ])
    zq = ArrayType(FT[
        5.72575636226887e-05,
        0.0001,
        5.72575636226887e-05,
        5.72575636226887e-05,
    ])
    z_rough = ArrayType(FT[Tuple(z0)[ii], Tuple(zt)[ii], Tuple(zq)[ii]])

    ts_sfc = TD.PhaseEquil_ρTq(param_set, ρ_sfc, T_sfc, qt_sfc)
    ts_in = TD.PhaseEquil_ρTq(param_set, ρ_ave, T_ave, qt_ave)

    u_in = u_ave

    # Constants
    Δz = Tuple(z)[jj]
    L_MO_init = FT(1.0)
    u_sfc = FT(0.0)
    args = (
        param_set,
        L_MO_init,
        ts_in, ts_sfc,
        u_in, u_sfc,
        z_rough,
        z_ave / 2,
    )

    #result = surface_conditions(args..., SF.DGScheme())
    result = surface_conditions(args..., SF.DGScheme())
    @show result
end
