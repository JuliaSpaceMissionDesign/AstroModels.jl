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
    compute_twobody(μ::T, pos::AbstractVector{T}) where T  

Compute acceleration due to the central body. Here `r` is the vector from the central body 
to the particle.
"""
@fastmath function compute_twobody(μ::Number, pos::AbstractVector{<:Number}) 
    r = SA[pos[1], pos[2], pos[3]]
    rn2 = r[1]*r[1] + r[2]*r[2] + r[3]*r[3]
    rn = sqrt(rn2)
    rn3 = rn2*rn 
    return -μ*r/rn3
end

"""
    compute_acceleration(mod::PointMass{T}, pos::AbstractVector{<:Number}, args...) where T 

Compute acceleration, given the position vector, `pos` of a spacecraft, due to a point mass.
"""
function compute_acceleration(mod::PointMass{T}, pos::AbstractVector{<:Number}, args...) where T 
    return compute_twobody(mod.μ, pos)
end