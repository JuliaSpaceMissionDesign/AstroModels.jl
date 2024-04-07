"""
    compute_thirdbody(R, Δi, μi) 

Compute acceleration due to the third body perturbation (only).
Here `R` is the vector from the central body to the particle, and `Δi` the
position vector from the central body to the perturbing body.

### References 
- Battin, R.H. -- An Introduction to the Mathematics and Methods of Astrodynamics, AIAA, 1999.
"""
@fastmath function compute_thirdbody(R, Δi, μi)

    @inbounds begin 
        Ri = R - Δi 

        # support variables 
        Δi2 = Δi[1]*Δi[1] + Δi[2]*Δi[2] + Δi[3]*Δi[3] 
        Δi3 = Δi2 * sqrt(Δi2)
        Ri2 = Ri[1]*Ri[1] + Ri[2]*Ri[2] + Ri[3]*Ri[3] 
        Ri3 = Ri2 * sqrt(Ri2)

        return - μi * Ri/Ri3 + μi * Δi/Δi3
    end

end

@fastmath function jacobian_thirdbody(R, Δi, μi)
    Ri = R - Δi
    return jacobian_twobody(μi, Ri)
end


"""
    compute_acceleration(center::PointMass{T}, pos::AbstractVector{N}, 
        third::AbstractVector{PointMass{T}}, axes, epoch, frames::G, 
        args...) where {T, N<:Number, G <:AbstractJSMDFrameGraph}

Compute accelerations due to third body perturbations in a planetocentric case.
"""
function compute_acceleration(center::PointMass{T}, pos::AbstractVector{N}, 
    third::AbstractVector{PointMass{T}}, axes, epoch, frames::G, 
    args...) where {T, N<:Number, G <:AbstractJSMDFrameGraph}

    a = compute_twobody(center.μ, pos)
    for p in third
        Δp = vector3(frames, center.id, p.id, axes, epoch)
        a += compute_thirdbody(pos, Δp, p.μ)
    end
    return a
    
end