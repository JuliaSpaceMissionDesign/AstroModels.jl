export AbstractAccelerationModel, AbstractAccelerationModelData, compute_acceleration

abstract type AbstractAccelerationModelData <: AbstractJSMDModelData end

abstract type AbstractAccelerationModel <: AbstractJSMDModel end

"""
    compute_acceleration(m::A, args...) where {A <: AbstractAccelerationModel}

This function serves as an interface for constructing acceleration models.

### Arguments
- `m::A`: An instance of a subtype of `AbstractAccelerationModel`, representing the 
    force model producing the acceleration.

!!! warning 
    Concrete implementations of `AbstractAccelerationModel` must provide this function!
"""
@interface function compute_acceleration(::A, args...) where {A <: AbstractAccelerationModel} end