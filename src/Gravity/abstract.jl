abstract type AbstractGravityModelData <: AbstractAccelerationModelData end

abstract type AbstractGravityHarmonicsData{N, T} <: AbstractGravityModelData end

abstract type AbstractGravityModel{T} <: AbstractAccelerationModel end