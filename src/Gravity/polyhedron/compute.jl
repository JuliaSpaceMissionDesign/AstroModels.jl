# TODO: implement 

@fastmath function compute_potential(p::PolyhedronGravity{T}, pvb, G, ρ, args...) where T
    r = SVector{3}(pvb[1], pvb[2], pvb[3])

    ue = zero(T)
    # Loop over edges 
    for e in p.edges    
        eᵢ, eⱼ = e.vertices 
        eᵢⱼ = e.magnitude 
        Eₑ = e.dyad
        @inbounds pᵢ, pⱼ = p.vertices[eᵢ], p.vertices[eⱼ]
        
        Rₑ, Rⱼ = pᵢ - r, pⱼ - r  
        rᵢ, rⱼ = norm(Rₑ), norm(Rⱼ)
        Lₑ = log((rᵢ + rⱼ + eᵢⱼ) / (rᵢ + rⱼ - eᵢⱼ))

        ue += Rₑ' * Eₑ * Rₑ * Lₑ
    end

    uf = zero(T)
    for f in p.faces 
        vᵢ, vⱼ, vₖ = f.vertices 
        @inbounds pᵢ, pⱼ, pₖ = p.vertices[vᵢ], p.vertices[vⱼ], p.vertices[vₖ]
        F = f.dyad

        Rᵢ, Rⱼ, Rₖ = pᵢ - r, pⱼ - r, pₖ - r   
        rᵢ, rⱼ, rₖ = norm(Rᵢ), norm(Rⱼ), norm(Rₖ)
        ωn = Rᵢ ⋅ (Rⱼ × Rₖ)
        ωd = rᵢ*rⱼ*rₖ + rᵢ*(Rⱼ⋅Rₖ) + rⱼ*(Rₖ⋅Rᵢ) + rₖ*(Rᵢ⋅Rⱼ)
        ω = 2*atan(ωn/ωd)

        uf += Rᵢ' * F * Rᵢ * ω
    end

    U = 0.5 * G * ρ * (ue - uf)
    return U
end

@fastmath function compute_acceleration(p::PolyhedronGravity{T}, pvb, G, ρ, args...) where T
    r = SVector{3}(pvb[1], pvb[2], pvb[3])

    δue = SVector{3, T}(0, 0, 0)
    # Loop over edges 
    for e in p.edges    
        eᵢ, eⱼ = e.vertices 
        eᵢⱼ = e.magnitude 
        Eₑ = e.dyad
        @inbounds pᵢ, pⱼ = p.vertices[eᵢ], p.vertices[eⱼ]
        
        Rₑ, Rⱼ = pᵢ - r, pⱼ - r  
        rᵢ, rⱼ = norm(Rₑ), norm(Rⱼ)
        Lₑ = log((rᵢ + rⱼ + eᵢⱼ) / (rᵢ + rⱼ - eᵢⱼ))

        δue += Eₑ * Rₑ * Lₑ
    end

    δuf = SVector{3, T}(0, 0, 0)
    for f in p.faces 
        vᵢ, vⱼ, vₖ = f.vertices 
        @inbounds pᵢ, pⱼ, pₖ = p.vertices[vᵢ], p.vertices[vⱼ], p.vertices[vₖ]
        F = f.dyad

        Rᵢ, Rⱼ, Rₖ = pᵢ - r, pⱼ - r, pₖ - r   
        rᵢ, rⱼ, rₖ = norm(Rᵢ), norm(Rⱼ), norm(Rₖ)
        ωn = Rᵢ ⋅ (Rⱼ × Rₖ)
        ωd = rᵢ*rⱼ*rₖ + rᵢ*(Rⱼ⋅Rₖ) + rⱼ*(Rₖ⋅Rᵢ) + rₖ*(Rᵢ⋅Rⱼ)
        ω = 2*atan(ωn/ωd)

        δuf += F * Rᵢ * ω
    end

    U = G * ρ * (- δue + δuf)
    return U
end

