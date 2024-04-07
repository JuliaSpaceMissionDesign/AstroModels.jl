using ForwardDiff
using LinearAlgebra
using StaticArrays
using Test

using AstroModels 
using AstroModels.Gravity


@testset "Jacobians" begin

    # point mass 
    μ = 1.0
    for i in 1:100
        r = randn(3) .+ 1

        @test_nowarn Gravity.compute_twobody(μ, r)
        @test_nowarn Gravity.jacobian_twobody(μ, r)

        j = Gravity.jacobian_twobody(μ, r)
        jr = ForwardDiff.jacobian(x->Gravity.compute_twobody(μ, x), r)
        @test all( isapprox(j, jr, atol=1e-13) )

    end

    # third body 
    μ = 1.0
    for i in 1:100
        v1 = randn(3) .+ 1 
        v2 = randn(3) .+ 1

        @test_nowarn Gravity.compute_thirdbody(v1, v2, μ)

        j = Gravity.jacobian_thirdbody(v1, v2, μ)
        jr = ForwardDiff.jacobian(x->Gravity.compute_thirdbody(x, v2, μ), v1)
        @test all( isapprox(j, jr, atol = 1e-13) )

    end
end;
