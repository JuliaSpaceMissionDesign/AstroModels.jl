using Test
using AstroModels 
using ForwardDiff

using LinearAlgebra
using StaticArrays
using Pkg.Artifacts

include("utils.jl");

@testset "Models" verbose=true begin
    @eval begin
        modules = [:Gravity]
        for m in modules 
            @testset "$m" verbose=true begin 
                include("$m/$m.jl")
            end 
        end
    end
end;
