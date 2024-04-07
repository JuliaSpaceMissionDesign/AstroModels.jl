export PointMass

"""
    PointMass{T} 

Gravitational model of a point mass.
"""
struct PointMass{T} <: AbstractGravityModel{T}
    id::Int
    μ::T 
end

"""
    parse_model(::Type{T}, ::Type{PointMass}, bodyid::Int, μ::T, args...) where {T}

Parse a [`PointMass`](@ref) model for a body with NAIFID `bodyid`, gravitational parameter `μ`
and `radius`.
"""
function parse_model(::Type{T}, ::Type{PointMass}, bodyid::Int, μ::T, args...) where {T}
    return PointMass{T}(bodyid, μ)
end


"""
    compute_twobody(μ::T, r::AbstractVector{T}) where T  

Compute acceleration due to the central body. Here `r` is the vector from the central body 
to the particle.
"""
@fastmath function compute_twobody(μ::Number, r::AbstractVector{<:Number}) 
    rn2 = r[1]*r[1] + r[2]*r[2] + r[3]*r[3]
    rn3 = rn2 * sqrt(rn2) 
    return -μ * r / rn3
end

@fastmath function jacobian_twobody(μ::Number, r::AbstractVector{T}) where T 
    
    @inbounds begin
        rn2 = r[1]*r[1] + r[2]*r[2] + r[3]*r[3]
        rn3 = rn2 * sqrt(rn2) 

        g = μ/rn3
        g2 = g/rn2

        Uxx = - g * fma(-1, 3*r[1]^2/rn2, 1)
        Uyy = - g * fma(-1, 3*r[2]^2/rn2, 1)
        Uzz = - g * fma(-1, 3*r[3]^2/rn2, 1)

        Uxy = g2 * 3*r[1]*r[2] 
        Uxz = g2 * 3*r[1]*r[3] 
        Uyz = g2 * 3*r[2]*r[3] 
    end

    return SMatrix{3, 3, T, 9}(
      Uxx,    Uxy,    Uxz, 
      Uxy,    Uyy,    Uyz,
      Uxz,    Uyz,    Uzz,     
    )'
end

"""
    compute_acceleration(mod::PointMass{T}, pos::AbstractVector{<:Number}, args...) where T 

Compute acceleration, given the position vector, `pos` of a spacecraft, due to a point mass.
"""
function compute_acceleration(mod::PointMass{T}, pos::AbstractVector{<:Number}, args...) where T 
    return compute_twobody(mod.μ, pos)
end