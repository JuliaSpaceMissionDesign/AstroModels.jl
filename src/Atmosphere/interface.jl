"""
    AbstractDensityModel{T}

Abstract type for all density models.
"""
abstract type AbstractDensityModel{T} <: AbstractJSMDModel end

"""
    compute()
"""
function compute_density end