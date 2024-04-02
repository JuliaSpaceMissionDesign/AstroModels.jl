
@inline @fastmath function fq(q)
    # Battin, page 389, eq. 8.60
    sq = sqrt(1+q)
    return q*((3 + q*(3 + q))/(1 + (sq*sq*sq)))
end

"""
    compute_thirdbody(μ::T, rc, rp) where T 

Compute acceleration due to the third body perturbation with an high-precision algorithm 
(no loss of significance results in the computations.)
Here `rc` is the vector from the central body to the particle, and `rp` the
position vector from the central body to the perturbing body.

### References 
- Battin, R.H. -- An Introduction to the Mathematics and Methods of Astrodynamics, AIAA, 1999.
"""
@fastmath function compute_thirdbody(μ::T, rc, rp) where T 

    r = SA[rc[1], rc[2], rc[3]]
    ρⱼ = SA[rp[1], rp[2], rp[3]]
    Δⱼ = r - ρⱼ
    
    r² = r[1]*r[1] + r[2]*r[2] + r[3]*r[3]
    ρ² = ρⱼ[1]*ρⱼ[1] + ρⱼ[2]*ρⱼ[2] + ρⱼ[3]*ρⱼ[3]
    Δ² = Δⱼ[1]*Δⱼ[1] + Δⱼ[2]*Δⱼ[2] + Δⱼ[3]*Δⱼ[3]
    Δ = sqrt(Δ²)
    Δ³ = Δ²*Δ

    # Battin, page 389, eq 8.62 
    δ = r[1]*ρⱼ[1] + r[2]*ρⱼ[2] + r[3]*ρⱼ[3]
    q = (r² - 2δ)/ρ²
    f = fq(q)

    # Battin, page 389, eq. 8.61
    return - μ/Δ³ * SA[
        f*ρⱼ[1] + r[1], 
        f*ρⱼ[2] + r[2],
        f*ρⱼ[3] + r[3]
    ]
end

"""
    compute_acceleration(center::PointMass{T}, pos::AbstractVector{<:Number}, third::PointMass{T},
        axes, epoch, frames::G, args...) where {T, G <:AbstractJSMDFrameGraph}

Compute accelerations due to third body perturbations given the position vector `pos` relative 
to the center, the current epoch `epoch` and a frame graph `frames`.
"""
function compute_acceleration(center::PointMass{T}, pos::AbstractVector{<:Number}, 
    third::PointMass{T}, axes, epoch, frames::G, args...) where {T, G <:AbstractJSMDFrameGraph}
    Δr = vector3(frames, center.id, third.id, axes, epoch)
    return compute_thirdbody(third.μ, pos, Δr)
end