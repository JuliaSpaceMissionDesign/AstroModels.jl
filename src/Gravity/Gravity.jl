module Gravity 

using Logging
using LinearAlgebra
using DelimitedFiles
using StaticArrays
using PreallocationTools

using JSMDInterfaces.Errors
using JSMDInterfaces.Frames: vector3, AbstractJSMDFrameGraph
using JSMDInterfaces.Interface: @interface
import JSMDInterfaces.Models: parse_data, parse_model

using JSMDUtils.Math: unitcross

using AstroModels: AbstractAccelerationModel, AbstractAccelerationModelData
import AstroModels: compute_acceleration

export parse_data, parse_model, compute_acceleration, compute_potential

# Interface types/methods
include("interface.jl")

# Point mass 
include("point/center.jl")
include("point/third.jl")

# Gravity harmonics
include("harmonics/parse.jl")
include("harmonics/type.jl")
include("harmonics/compute.jl")

# Polyhedron 
include("polyhedron/parse.jl")
include("polyhedron/type.jl")
include("polyhedron/compute.jl")

end