export PointMass, parse_model, compute_acceleration 

struct PointMass{T} <: AbstractGravityModel{T}
    bodyid::Int
    μ::T # L³/T²
    radius::T # L
end

"""
    parse_model(::Type{T}, ::Type{PointMass}, bodyid::Int, μ::T, radius::T, args...) where {T}

Parse a [`PointMass`](@ref) model for a body with NAIFID `bodyid`, gravitational parameter `μ`
and `radius`.
"""
function parse_model(::Type{T}, ::Type{PointMass}, bodyid::Int, μ::T, radius::T, args...) where {T}
    return PointMass{T}(bodyid, μ, radius)
end


"""
    compute_twobody(μ::T, pos::AbstractVector{T}) where T  

Compute acceleration due to the central body. Here `r` is the vector from the central body 
to the particle.
"""
@fastmath function compute_twobody(μ::T, pos::AbstractVector{T}) where T 
    r = SA[pos[1], pos[2], pos[2]]
    rn2 = r[1]*r[1] + r[2]*r[2] + r[3]*r[3]
    rn = sqrt(rn2)
    rn3 = rn2*rn 
    return -μ*r/rn3
end

function compute_acceleration(mod::PointMass{T}, pos::AbstractVector{T}, args...) where T 
    return compute_twobody(mod.μ, pos)
end