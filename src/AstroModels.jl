module AstroModels

    using JSMDInterfaces.Models: AbstractJSMDModelData, AbstractJSMDModel
    using JSMDInterfaces.Interface: @interface
    using PreallocationTools
    
    include("interface.jl")
    include(joinpath("Gravity", "Gravity.jl"))
    include(joinpath("SolarPressure", "SolarPressure.jl"))
    include(joinpath("Atmosphere", "Atmosphere.jl"))

end