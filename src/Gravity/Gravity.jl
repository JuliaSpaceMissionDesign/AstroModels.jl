module Gravity 

using Logging
using DelimitedFiles
using StaticArrays
using PreallocationTools

using JSMDInterfaces.Errors
using JSMDInterfaces.Frames: vector3, AbstractJSMDFrameGraph
import JSMDInterfaces.Models: parse_data, parse_model

using AstroModels: AbstractAccelerationModel, AbstractAccelerationModelData

# Abstract types/methods
include("abstract.jl")

# Gravity harmonics
include("harmonics/parse.jl")
include("harmonics/type.jl")
include("harmonics/compute.jl")

# Point mass 
include("point/center.jl")
include("point/third.jl")

end