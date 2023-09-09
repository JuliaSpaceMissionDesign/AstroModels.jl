
"""
    compute_thirdbody(μ::T, r::AbstractVector{T}, d::AbstractVector{T}) where T 

Compute acceleration due to the third body perturbation.
Here `r` is the vector from the central body to the particle, and `d`  the
position vector from the central body to the perturbing body.
"""
@fastmath function compute_thirdbody(μ::T, r::AbstractVector{T}, d::AbstractVector{T}) where T 
    @inbounds @views begin
        D = d[1:3]
        Δ = r[1:3] - D
        Dn = sqrt(D[1]*D[1] + D[2]*d[2] + D[3]*D[3])
        Δn = sqrt(Δ[1]*Δ[1] + Δ[2]*Δ[2] + Δ[3]*Δ[3])
    end
    return - μ * (Δ/Δn^3 + D/Dn^3)
end