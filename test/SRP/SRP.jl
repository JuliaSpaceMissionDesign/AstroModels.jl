using ForwardDiff
using LinearAlgebra
using StaticArrays
using Test

using AstroModels 
using AstroModels.SolarPressure

@testset "Models consistency" begin
    ρs, ρd = 1.0, 0.0
    A, m, P = 1.0, 1.0, 1.0

    for _ in 1:1000
        v = @SVector randn(3)
        v /= norm(v)

        ac = compute_srp_cannonball(ρs, A, m, P, v)
        af = compute_srp_flatplate(ρs, ρd, A, m, P, -v, v)
        @test maximum(abs.(ac - af)) ≈ 0.0 atol=1e-14
    end
end

@testset "Flatplate" begin
    ρs, ρd = 1.0, 0.0
    A, m, P = 1.0, 1.0, 1.0

    # Plate inclined by 60 degrees 
    vp = [cosd(60), sind(60), 0]
    v = [1., 0., 0.]

    ac = compute_srp_cannonball(ρs, A, m, P, v)
    af = compute_srp_flatplate(ρs, ρd, A, m, P, -vp, v)

    @test ac[1]/8 ≈ af[1] atol=1e-15
    @test af[2]/sqrt(3) ≈ af[1]
    @test ac[3] ≈ 0.0

    # Plate at 90 degrees 
    af = compute_srp_flatplate(ρs, ρd, A, m, P, [0., 1., 0], v)
    @test all( af .≈ 0 )

    # Plate at ϕ
    ρs, ρd = 0.5, 0.1
    for ϕ in -π/3:π/10:π/3
        sϕ, cϕ = sincos(ϕ)
        vp = [-cϕ, sϕ, 0]
        cθ = dot(vp, v)

        af = compute_srp_flatplate(ρs, ρd, A, m, P, vp , v)
        ρs, ρd = 0.5, 0.1

        f = (1-ρs - 2*cϕ*(ρs*cθ + ρd/3)) / (2*sϕ*(ρs*cθ + ρd/3))
        
        @test af[2] * f ≈ af[1] atol=1e-14
        @test af[3] == 0
    end
    
end;

@testset "Shadow models" verbose=true begin
    include("shadow.jl")
end