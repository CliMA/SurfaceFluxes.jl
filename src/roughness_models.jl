abstract type RoughnessModel end
struct CharnockRoughness <: RoughnessModel end
struct ScalarRoughness <: RoughnessModel end

# Roughness Model Computations
function compute_z0(u★, sfc_param_set, 
                    sc::AbstractSurfaceConditions, ::ScalarRoughness, ::UF.MomentumTransport)
    return sc.z0m
end

function compute_z0(u★, sfc_param_set, 
                    sc, ::CharnockRoughness, ::UF.MomentumTransport)
    𝛼 = eltype(u★)(0.011)
    𝑔 = SFP.grav(sfc_param_set)
    return 𝛼 * u★^2 / 𝑔
end
function compute_z0(u★, sfc_param_set, 
                    sc, ::CharnockRoughness, ::UF.HeatTransport)
    𝛼 = eltype(u★)(0.011)
    𝑔 = SFP.grav(sfc_param_set)
    return 𝛼 * u★^2 / 𝑔
end

function compute_z0(u★, sfc_param_set, 
                    sc::AbstractSurfaceConditions, ::ScalarRoughness, ::UF.HeatTransport)
    return sc.z0m
end
function compute_z0(u★, sfc_param_set, 
                    sc::Coefficients, ::ScalarRoughness, ::UF.MomentumTransport)
    return compute_z0(u★, sfc_param_set, sc, CharnockRoughness(), UF.MomentumTransport())
end
function compute_z0(u★, sfc_param_set, 
                    sc::Coefficients, ::ScalarRoughness, ::UF.HeatTransport)
    return compute_z0(u★, sfc_param_set, sc, CharnockRoughness(), UF.MomentumTransport())
end
function compute_z0(u★, sfc_param_set, 
                    sc, ::ScalarRoughness, ::UF.HeatTransport)
    return sc.z0b
end

