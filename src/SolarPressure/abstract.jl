abstract type AbstractSolarPressureModelData{T} <: AbstractAccelerationModelData end

abstract type AbstractSolarPressureModel{T} <: AbstractAccelerationModel end

struct NoSolarPressureModel{T} <: AbstractSolarPressureModel{T} end 

@fastmath function compute_acceleration(::NoSolarPressureModel{T}, args...; kwargs...) where T 
    return SA[0., 0., 0.]
end