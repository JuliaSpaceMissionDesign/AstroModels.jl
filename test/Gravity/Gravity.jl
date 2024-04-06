using ForwardDiff
using LinearAlgebra
using StaticArrays
using Pkg.Artifacts
using Test

using AstroModels 
using AstroModels.Gravity

@time @testset "Spherical Harmonics" verbose=true begin
    include("harmonics.jl")
end;

@time @testset "Polyhedron" verbose=true begin 
    include("polyhedron.jl")
end;
