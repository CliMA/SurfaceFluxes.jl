using NCDatasets
using Plots
using LaTeXStrings

data = NCDataset("DYCOMS_RF01.nc");
z = data["z_half"];
ρ = data.group["profiles"]["rho"];
u = data.group["profiles"]["u_mean"];
v = data.group["profiles"]["v_mean"];
qt = data.group["profiles"]["qt_mean"];
ql = data.group["profiles"]["ql_mean"];
p0 = data.group["profiles"]["p0"];
b = data.group["profiles"]["buoyancy_mean"];
θli = data.group["profiles"]["thetali_mean"];

function mean(X)
    return mean(X, dims=2)
end
