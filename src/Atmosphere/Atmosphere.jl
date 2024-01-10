module Atmosphere 

using AstroModels: AbstractAccelerationModel, AbstractAccelerationModelData
import JSMDInterfaces.Models: parse_data, parse_model
using JSMDInterfaces.Errors
using StaticArrays

include("abstract.jl")
include("exp.jl")

end