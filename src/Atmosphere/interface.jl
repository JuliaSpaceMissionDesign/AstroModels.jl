"""
    AbstractDensityModel{T}

Abstract type for all density models.
"""
abstract type AbstractDensityModel{T} <: AbstractJSMDModel end

"""
    compute_density(model, pos)

Compute density based on a given model in a given position.
"""
function compute_density end

"""
    compute_wind(model, pos)

Compute wind based on a given model in a given position.
"""
function compute_wind end