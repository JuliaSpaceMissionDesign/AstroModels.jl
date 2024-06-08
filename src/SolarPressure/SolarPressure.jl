module SolarPressure

using StaticArrays
using LinearAlgebra

using JSMDInterfaces.Models: AbstractJSMDModel
using JSMDUtils.Math: unitvec
using AstroModels: AbstractAccelerationModel, AbstractAccelerationModelData

include("interface.jl")
include("pressure.jl")
include("cannonball.jl")
include("flatplate.jl")
include("shadow.jl")

end