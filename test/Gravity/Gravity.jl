using ForwardDiff
using LinearAlgebra
using StaticArrays
using Pkg.Artifacts
using Test

using AstroModels 
using AstroModels.Gravity

@testset "Spherical Harmonics" verbose=true begin
    include("harmonics.jl")
end;

@testset "Polyhedron" verbose=true begin 
    include("polyhedron.jl")
end;
