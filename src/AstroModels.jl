module AstroModels

    using JSMDInterfaces.Models: AbstractJSMDModelData, AbstractJSMDModel
    using JSMDInterfaces.Interface: @interface
    using PreallocationTools
    
    include("abstract.jl")

    include("Gravity/Gravity.jl") 
    include("SolarPressure/SolarPressure.jl") 
    include("Atmosphere/Atmosphere.jl") 
end