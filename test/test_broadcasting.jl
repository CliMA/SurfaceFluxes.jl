
using Test
import SurfaceFluxes as SF

# Check that broadcasting works as intended (treating params as scalars)
@testset "Parameter Broadcasting" begin
    FT = Float64

    # define dummy inputs
    N = 10
    x = rand(FT, N)

    # 1. Roughness Params
    roughness = SF.ConstantRoughnessParams(z0m = FT(0.001), z0s = FT(0.0001))

    # Dummy function that mimics usage
    f(val, p::SF.AbstractRoughnessParams) = val * p.z0m

    # Broadcast
    res = f.(x, roughness)

    @test length(res) == N
    @test res ≈ x .* FT(0.001)

    # 2. Gustiness Specs
    gustiness = SF.ConstantGustinessSpec(FT(1.0))
    g(val, p::SF.AbstractGustinessSpec) = val * p.value

    res_g = g.(x, gustiness)

    @test length(res_g) == N
    @test res_g ≈ x .* FT(1.0)

    # 3. Check COARE3
    coare = SF.COARE3RoughnessParams{FT}()
    h(val, p::SF.COARE3RoughnessParams) = val * p.α_low

    res_h = h.(x, coare)
    @test length(res_h) == N
    @test res_h ≈ x .* coare.α_low
end
