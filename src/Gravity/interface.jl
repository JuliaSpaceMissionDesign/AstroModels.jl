export compute_potential

"""
    AbstractGravityModel

Abstract type for all data associated to gravitational models.
"""
abstract type AbstractGravityModelData <: AbstractAccelerationModelData end

"""
    AbstractGravityModel{T}

Abstract type for all gravitational models.
"""
abstract type AbstractGravityModel{T} <: AbstractAccelerationModel end

"""
    compute_potential(m::A, args...) where {A <: AbstractAccelerationModel}

This function serves as an interface for constructing gravitational potentials models.

### Arguments
- `m`: An instance of a subtype of `AbstractGravityModel`, representing the 
    gravity model producing the potential.

!!! warning 
    Concrete implementations of `AbstractGravityModel` must provide this function!
"""
@interface function compute_potential(::A, args...) where {A <: AbstractGravityModel} end
