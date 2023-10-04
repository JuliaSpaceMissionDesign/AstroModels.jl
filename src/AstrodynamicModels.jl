module AstrodynamicModels

    using JSMDInterfaces.Models: AbstractJSMDModelData, AbstractJSMDModel
    include("abstract.jl")

    include("Gravity/Gravity.jl") 
    include("SolarPressure/SolarPressure.jl") 
    include("Atmosphere/Atmosphere.jl") 
end