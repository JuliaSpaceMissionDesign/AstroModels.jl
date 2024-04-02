export AbstractAccelerationModel, AbstractAccelerationModelData, compute_acceleration

abstract type ForceModelData <: AbstractJSMDModelData end

abstract type ForceModel <: AbstractJSMDModel end

abstract type AbstractAccelerationModelData <: ForceModelData end

abstract type AbstractAccelerationModel <: ForceModel end

"""
    compute_acceleration(m::A, args...) where {A <: AccelerationModel}

This function serves as an interface for constructing acceleration models.

### Arguments
- `m::A`: An instance of a subtype of `AccelerationModel`, representing the force model 
    producing the acceleration.

!!! warning 
    Concrete implementations of `AbstractAccelerationModel` must provide this function!
"""
@interface function compute_acceleration(::A, args...) where {A <: AbstractAccelerationModelData} end