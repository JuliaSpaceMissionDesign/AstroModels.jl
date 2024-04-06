module SolarPressure

using StaticArrays
using JSMDInterfaces.Models: AbstractJSMDModel
using AstroModels: AbstractAccelerationModel, AbstractAccelerationModelData

include("interface.jl")
include("pressure.jl")
include("cannonball.jl")
include("flatplate.jl")

end