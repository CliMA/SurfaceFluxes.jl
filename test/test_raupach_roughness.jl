# Tests for Raupach (1994) canopy roughness parameterization

module TestRaupachRoughness

using Test
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions as UF
import ClimaParams as CP

@testset "Raupach Roughness Parameterization" begin
    FT = Float64
    param_set = SFP.SurfaceFluxesParameters(FT, UF.BusingerParams)

    # Default Raupach parameters
    spec = SF.RaupachRoughnessParams{FT}()
    u_star = FT(0.3)
    # Roughness-sublayer influence function at the canopy top (Raupach 1994, Eq. 5)
    Ψ_h = log(spec.c_w) - 1 + 1 / spec.c_w

    @testset "ClimaParams defaults" begin
        # The TOML defaults are those of the struct
        @test SF.RaupachRoughnessParams(CP.create_toml_dict(FT)) == spec
        @test Ψ_h ≈ 0.193 atol = 1e-3
    end

    @testset "Basic Functionality" begin
        # Normal canopy: PAI = 3, h = 10m
        inputs = (PAI = FT(3.0), h = FT(10.0))

        z0m = SF.momentum_roughness(spec, u_star, param_set, inputs)
        z0s = SF.scalar_roughness(spec, u_star, param_set, inputs)

        @test z0m > 0
        @test z0s > 0
        @test z0s < z0m  # Scalar roughness < momentum roughness (Stanton number < 1)
        @test z0m < inputs.h  # Roughness length < canopy height
    end

    @testset "PAI Dependence" begin
        h = FT(10.0)

        # Raupach formula is non-monotonic in PAI: z0m increases with PAI at low
        # PAI (λ < 0.29, PAI < 0.58) and decreases at higher PAI once u_star/U(h)
        # reaches its sheltering cap of 0.3 (Eq. 8 in Raupach 1994)
        inputs_sparse = (PAI = FT(0.1), h = h)
        inputs_peak = (PAI = FT(0.6), h = h)
        inputs_mid = (PAI = FT(2.0), h = h)
        inputs_dense = (PAI = FT(4.0), h = h)

        z0m_sparse = SF.momentum_roughness(spec, u_star, param_set, inputs_sparse)
        z0m_peak = SF.momentum_roughness(spec, u_star, param_set, inputs_peak)
        z0m_mid = SF.momentum_roughness(spec, u_star, param_set, inputs_mid)
        z0m_dense = SF.momentum_roughness(spec, u_star, param_set, inputs_dense)

        @test z0m_peak > z0m_sparse
        @test z0m_peak > z0m_mid > z0m_dense

        # All values should be reasonable fractions of canopy height
        @test z0m_sparse < h
        @test z0m_peak < h
    end


    @testset "Canopy Height Dependence" begin
        PAI = FT(3.0)

        # Taller canopy => Higher z0m
        inputs_short = (PAI = PAI, h = FT(5.0))
        inputs_tall = (PAI = PAI, h = FT(20.0))

        z0m_short = SF.momentum_roughness(spec, u_star, param_set, inputs_short)
        z0m_tall = SF.momentum_roughness(spec, u_star, param_set, inputs_tall)

        @test z0m_tall > z0m_short
    end

    @testset "Floor on the frontal area index" begin
        # With a floor, at PAI = 0 (for example, a leaf area index without stems) the
        # frontal area index is λ_min, so a tall canopy keeps a roughness length far
        # above the fixed value
        spec_floor = SF.RaupachRoughnessParams{FT}(λ_min = 0.05)
        inputs_zero = (PAI = FT(0.0), h = FT(10.0))
        inputs_floor =
            (PAI = spec_floor.λ_min / spec_floor.frontal_area_ratio, h = FT(10.0))
        z0m_zero = SF.momentum_roughness(spec_floor, u_star, param_set, inputs_zero)
        z0m_floor = SF.momentum_roughness(spec_floor, u_star, param_set, inputs_floor)

        @test z0m_zero ≈ z0m_floor
        @test z0m_zero > 100 * SFP.z0m_fixed(param_set)
        @test SF.displacement_height(spec_floor, inputs_zero) ≈
              SF.displacement_height(spec_floor, inputs_floor)

        # Without a canopy, the roughness length is the fixed value
        inputs_bare = (PAI = FT(0.0), h = FT(0.0))
        @test SF.momentum_roughness(spec, u_star, param_set, inputs_bare) ≈
              SFP.z0m_fixed(param_set)
        @test SF.displacement_height(spec, inputs_bare) == 0

        # By default there is no floor: a canopy with zero plant area index has no
        # displacement, and its roughness length is that of the substrate drag alone
        @test spec.λ_min == 0
        @test SF.raupach_displacement_fraction(spec, FT(0)) ≈ 0 atol = 1e-6
        @test SF.displacement_height(spec, inputs_zero) ≈ 0 atol = 1e-5
        κ = SFP.von_karman_const(param_set)
        @test SF.momentum_roughness(spec, u_star, param_set, inputs_zero) ≈
              max(
            inputs_zero.h * exp(-κ / sqrt(spec.C_S) - Ψ_h),
            SFP.z0m_fixed(param_set),
        ) rtol = 1e-6
    end

    @testset "Closed-form values" begin
        # λ = 0.5 (plant area index Λ = 1 at frontal_area_ratio = 0.5), above the
        # sheltering cap; d / h depends on the canopy area index Λ (Raupach 1994,
        # Eq. 8, and Fig. 1b, where d / h ≈ 0.66 at Λ = 1)
        κ = SFP.von_karman_const(param_set)
        h = FT(10)
        inputs = (PAI = FT(1.0), h = h)
        λ = FT(0.5)
        Λ = FT(1.0)
        @test SF.frontal_area_index(spec, inputs.PAI) ≈ λ
        @test SF.canopy_area_index(spec, inputs.PAI) ≈ Λ
        x = sqrt(FT(7.5) * Λ)
        d_over_h = 1 - (1 - exp(-x)) / x
        z0m_over_h = (1 - d_over_h) * exp(-κ / FT(0.3) - Ψ_h)
        @test SF.displacement_height(spec, inputs) ≈ d_over_h * h
        @test SF.momentum_roughness(spec, u_star, param_set, inputs) ≈
              z0m_over_h * h
        @test d_over_h ≈ 0.6585 atol = 1e-3
        @test z0m_over_h ≈ 0.0742 atol = 1e-3

        # Below the cap, u★ / U(h) = sqrt(C_S + C_R λ)
        inputs_sparse = (PAI = FT(0.2), h = h)
        λ_sparse = FT(0.1)
        x_sparse = sqrt(FT(7.5) * inputs_sparse.PAI)
        d_sparse = 1 - (1 - exp(-x_sparse)) / x_sparse
        z0m_sparse =
            (1 - d_sparse) *
            exp(-κ / sqrt(FT(0.003) + FT(0.3) * λ_sparse) - Ψ_h)
        @test SF.momentum_roughness(spec, u_star, param_set, inputs_sparse) ≈
              z0m_sparse * h
    end

    @testset "Displacement height" begin
        h = FT(20)
        d = [
            SF.displacement_height(spec, (PAI = PAI, h = h)) for
            PAI in FT.((0.5, 1, 2, 4, 8))
        ]
        # d / h rises monotonically with PAI toward the canopy top
        @test all(diff(d) .> 0)
        @test all(0 .< d .< h)
        @test d[end] / h > 0.8
        # d + z0m stays below the canopy top
        for PAI in FT.((0.1, 0.6, 2, 8))
            inputs = (PAI = PAI, h = h)
            z0m = SF.momentum_roughness(spec, u_star, param_set, inputs)
            @test SF.displacement_height(spec, inputs) + z0m < h
        end
    end

    @testset "Frontal area ratio" begin
        inputs = (PAI = FT(2.0), h = FT(10.0))
        spec_half = SF.RaupachRoughnessParams{FT}(frontal_area_ratio = 0.25)
        inputs_half = (PAI = FT(1.0), h = FT(10.0))
        # The ratio enters the drag partition only: halving it halves the frontal
        # area index, while the displacement height follows the plant area index
        @test SF.frontal_area_index(spec_half, inputs.PAI) ≈
              SF.frontal_area_index(spec, inputs_half.PAI)
        @test SF.displacement_height(spec_half, inputs) ≈
              SF.displacement_height(spec, inputs)
        κ = SFP.von_karman_const(param_set)
        @test SF.canopy_area_index(spec_half, inputs.PAI) ≈ inputs.PAI
        @test SF.momentum_roughness(spec_half, u_star, param_set, inputs) ≈
              inputs.h * SF.raupach_roughness_fraction(spec_half, κ, inputs.PAI)
        # With the floor active, both indices describe the same canopy
        spec_floor = SF.RaupachRoughnessParams{FT}(λ_min = 0.05)
        inputs_zero = (PAI = FT(0.0), h = FT(10.0))
        @test SF.canopy_area_index(spec_floor, inputs_zero.PAI) ≈
              spec_floor.λ_min / spec_floor.frontal_area_ratio
        # Integer plant area indices give the same results as floating-point ones
        @test SF.frontal_area_index(spec, 2) === SF.frontal_area_index(spec, 2.0)
        @test SF.canopy_area_index(spec, 2) ≈ 2
        @test SF.displacement_height(spec, (PAI = 2, h = 10)) ≈
              SF.displacement_height(spec, (PAI = 2.0, h = 10.0))
        @test SF.momentum_roughness(spec, u_star, param_set, (PAI = 2, h = 10)) ≈
              SF.momentum_roughness(spec, u_star, param_set, (PAI = 2.0, h = 10.0))
        @test SF.raupach_displacement_fraction(spec, 1) ≈
              SF.raupach_displacement_fraction(spec, 1.0)
        # The deprecated field name `LAI` gives the same results as `PAI`
        @test SF.displacement_height(spec, (LAI = 2.0, h = 10.0)) ===
              SF.displacement_height(spec, (PAI = 2.0, h = 10.0))
        @test SF.momentum_roughness(spec, u_star, param_set, (LAI = 2.0, h = 10.0)) ===
              SF.momentum_roughness(spec, u_star, param_set, (PAI = 2.0, h = 10.0))
        # Float64 coefficients keep Float32 inputs in Float32
        spec64 = SF.RaupachRoughnessParams{Float64}()
        inputs32 = (PAI = 1.0f0, h = 10.0f0)
        @test SF.displacement_height(spec64, inputs32) isa Float32
        @test SF.raupach_roughness_fraction(spec64, 0.4f0, 1.0f0) isa Float32
    end

    @testset "Stanton Number Scaling" begin
        inputs = (PAI = FT(3.0), h = FT(10.0))

        z0m, z0s = SF.momentum_and_scalar_roughness(spec, u_star, param_set, inputs)

        # z0s = z0m * stanton_number
        @test z0s ≈ z0m * spec.stanton_number

        # kB⁻¹ = 2 via the Stanton number
        spec_kB = SF.RaupachRoughnessParams{FT}(stanton_number = exp(-FT(2)))
        z0m_kB, z0s_kB =
            SF.momentum_and_scalar_roughness(spec_kB, u_star, param_set, inputs)
        @test z0m_kB ≈ z0m
        @test log(z0m_kB / z0s_kB) ≈ 2
    end

    @testset "Float32" begin
        spec32 = SF.RaupachRoughnessParams{Float32}()
        param_set32 = SFP.SurfaceFluxesParameters(Float32, UF.BusingerParams)
        inputs32 = (PAI = 2.0f0, h = 10.0f0)
        @test SF.momentum_roughness(spec32, 0.3f0, param_set32, inputs32) isa Float32
        @test SF.displacement_height(spec32, inputs32) isa Float32
    end
end

end # module
