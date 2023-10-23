export compute_acceleration

using FrameTransformations.Frames: vector3, FrameSystem

struct ThirdBodies{T} <: AbstractGravityModel{T}
    cid::Int
    axes::Int
    vid::Vector{Int}
    μs::Vector{T}
end

"""
    parse_model(::Type{T}, ::Type{ThirdBodies}, bodyid::Int, μ::T, radius::T, args...) where {T}

Parse a [`ThirdBodies`](@ref) model representing third body perturbations for the bodies with 
NAIFIDs `vid` and gravitational parameter `μs` when computing the acceleration from relative 
to a central body `cid` w.r.t. the `axes` frame.
"""
function parse_model(::Type{T}, ::Type{ThirdBodies}, cid::Int, axes::Int, vid::Vector{Int}, μs::Vector{T}, args...) where {T}
    return ThirdBodies{T}(cid, axes, vid, μs)
end


"""
    compute_thirdbody(μ::T, r::AbstractVector{T}, d::AbstractVector{T}) where T 

Compute acceleration due to the third body perturbation.
Here `r` is the vector from the central body to the particle, and `p` the
position vector from the central body to the perturbing body.
"""
@fastmath function compute_thirdbody(μ::T, r::AbstractVector{T}, p::AbstractVector{T}) where T 
    ri = SA[p[1], p[2], p[3]]
    Δ = SA[r[1], r[2], r[3]] - ri
    rin = sqrt(ri[1]*ri[1] + ri[2]*ri[2] + ri[3]*ri[3])
    Δn = sqrt(Δ[1]*Δ[1] + Δ[2]*Δ[2] + Δ[3]*Δ[3])
    rin3 = rin*rin*rin
    Δn3 = Δn*Δn*Δn 
    return -μ * SA[
        ri[1]/rin3 + Δ[1]/Δn3, 
        ri[2]/rin3 + Δ[2]/Δn3, 
        ri[3]/rin3 + Δ[3]/Δn3
    ]
end

@inline @fastmath function fq(q)
    # Battin, page 389, eq. 8.60
    sq = sqrt(1+q)
    return q*((3 + q*(3 + q))/(1 + (sq*sq*sq)))
end

"""
    compute_thirdbody_hp(μ::T, r::AbstractVector{T}, d::AbstractVector{T}) where T 

Compute acceleration due to the third body perturbation with an high-precision algorithm 
(no loss of significance results in the computations.)
Here `rc` is the vector from the central body to the particle, and `rp` the
position vector from the central body to the perturbing body.

### References 
- Battin, R.H. -- An Introduction to the Mathematics and Methods of Astrodynamics, AIAA, 1999.
"""
@fastmath function compute_thirdbody_hp(μ::T, rc::AbstractVector{T}, rp::AbstractVector{T}) where T 

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

function compute_acceleration(mod::ThirdBodies{T}, pos::AbstractVector{T}, epoch, 
    frames::FrameSystem{N, T}, args...) where {N, T} 
    acc_tot = SVector{3, T}(0., 0., 0.)
    # Loop over the perturbing bodies
    for (pid, μi) in zip(mod.vid, mod.μs)
        Δr = vector3(frames, mod.cid, pid, mod.axes, epoch)
        acc_tot += compute_thirdbody_hp(μi, pos, Δr)
    end
    return acc_tot
end