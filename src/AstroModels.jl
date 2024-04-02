module AstroModels

    using JSMDInterfaces.Models: AbstractJSMDModelData, AbstractJSMDModel
    using JSMDInterfaces.Interface: @interface
    using PreallocationTools
    
    include("interface.jl")
    include("Gravity/Gravity.jl") 

end