abstract type AbstractGravityModelData <: AbstractJSMDModelData end

abstract type AbstractGravityHarmonicsData{N, T} <: AbstractGravityModelData end

abstract type AbstractGravityModel{T} <: AbstractJSMDModel end