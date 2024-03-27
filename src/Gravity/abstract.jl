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
    AbstractGravityHarmonicsData{N, T}

Abstract type for all spherical harmonics expansion model data.
"""
abstract type AbstractGravityHarmonicsData{N, T} <: AbstractGravityModelData end
