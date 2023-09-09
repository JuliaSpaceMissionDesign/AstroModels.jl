@fastmath function compute_cannonball_srp(ρ::T, P::T, Asc::T, Msc::T, 
    r::AbstractVector{T}, s::AbstractVector{T}) where T 
    @inbounds @views begin
        Δ = s[1:3] - r[1:3]
        Δn = sqrt(Δ[1]*Δ[1] + Δ[2]*Δ[2] + Δ[3]*Δ[3])
    end
    return - (1+ρ)*P*Asc/Msc * Δ/Δn^3
end