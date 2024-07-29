
struct CIRA72{T} <: AbstractDensityModel{T} 
    h::Vector{T}
    ρ::Vector{T}
    S::Vector{T}
    T::Vector{T}
    w::Vector{T}
end

function CIRA72{T}() where T
    return CIRA72{T}(
        convert.(T, _CIRA72_h), convert.(T, _CIRA72_d), 
        convert.(T, _CIRA72_S), convert.(T, _CIRA72_T), 
        convert.(T, _CIRA72_w)
    )
end

function interpolate(value1, value2, fraction)
    return value1 + (value2 - value1) * fraction
end

function compute_density(m::CIRA72{<:Number}, h::Number) # h in km, ρ in kg/m³
    idx = findlast(x -> x <= h, m.h)
    Δh = h - m.h[idx]
    return m.ρ[idx] * exp( - Δh/m.S[idx] )
end