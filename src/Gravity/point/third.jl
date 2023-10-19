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
    ri = SA[p[1], p[2], p[2]]
    Δ = SA[r[1], r[2], r[2]] - ri
    rin = sqrt(ri[1]*ri[1] + ri[2]*ri[2] + ri[3]*ri[3])
    Δn = sqrt(Δ[1]*Δ[1] + Δ[2]*Δ[2] + Δ[3]*Δ[3])
    rin3 = rin*rin*rin
    Δn3 = Δn*Δn*Δn 
    return -μ * SVector{3, T}(
        ri[1]/rin3 + Δ[1]/Δn3, 
        ri[2]/rin3 + Δ[2]/Δn3, 
        ri[3]/rin3 + Δ[3]/Δn3
    )
end

"""
    compute_thirdbody_hp(μ::T, r::AbstractVector{T}, d::AbstractVector{T}) where T 

Compute acceleration due to the third body perturbation with an high-precision algorithm 
(no loss of significance results in the computations.)
Here `r` is the vector from the central body to the particle, and `p` the
position vector from the central body to the perturbing body.

### References 
- Battin, R.H. -- An Introduction to the Mathematics and Methods of Astrodynamics, AIAA, 1999.
"""
@fastmath function compute_thirdbody_hp(μ::T, r::AbstractVector{T}, p::AbstractVector{T}) where T 
    
    ri = SA[p[1], p[2], p[2]]
    ri2 = ri[1]*ri[1] + ri[2]*ri[2] + ri[3]*ri[3]
    
    Δ = SA[r[1], r[2], r[2]] - ri
    Δn2 = Δ[1]*Δ[1] + Δ[2]*Δ[2] + Δ[3]*Δ[3]
    Δn = sqrt(Δn2)
    Δn3 = Δn2 * Δn
    
    # Battin, page 389, eq 8.59/8.60
    Δm = -(Δ + ri)
    q = (r[1]*Δm[1] + r[2]*Δm[2] + r[3]*Δm[3])/ri2
    tmp = sqrt(1+q)
    f = tmp*(1+q) - 1 # f(q)

    # Battin, page 389, eq. 8.61
    return μ/Δn3 * SVector{3, T}(
        f*ri[1] + r[1], 
        f*ri[2] + r[2],
        f*ri[3] + r[3]
    )
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