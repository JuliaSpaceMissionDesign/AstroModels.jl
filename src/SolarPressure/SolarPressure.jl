module SolarPressure

using StaticArrays

using JSMDInterfaces.Models: AbstractJSMDModel
using AstroModels: AbstractAccelerationModel, AbstractAccelerationModelData
import JSMDInterfaces.Models: parse_data, parse_model
using JSMDUtils.Math: unitvec

include("abstract.jl")
include("pressure.jl")
include("cannonball.jl")
include("flatplate.jl")

end