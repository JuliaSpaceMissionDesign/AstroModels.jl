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
