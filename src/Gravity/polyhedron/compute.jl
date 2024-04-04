
@fastmath function precompute_edges(model, edge, r)
    eᵢ, eⱼ = edge.vertices 
    eᵢⱼ = edge.magnitude 
    Eₑ = edge.dyad
    @inbounds pᵢ, pⱼ = model.vertices[eᵢ], model.vertices[eⱼ]
    
    Rₑ, Rⱼ = pᵢ - r, pⱼ - r  
    rᵢ, rⱼ = norm(Rₑ), norm(Rⱼ)
    @fastmath Lₑ = log((rᵢ + rⱼ + eᵢⱼ) / (rᵢ + rⱼ - eᵢⱼ))
    return Rₑ, Eₑ, Lₑ
end

@fastmath function precompute_faces(model, face, r)
    vᵢ, vⱼ, vₖ = face.vertices 
    F = face.dyad
    @inbounds pᵢ, pⱼ, pₖ = model.vertices[vᵢ], model.vertices[vⱼ], model.vertices[vₖ]

    Rᵢ, Rⱼ, Rₖ = pᵢ - r, pⱼ - r, pₖ - r   
    rᵢ, rⱼ, rₖ = norm(Rᵢ), norm(Rⱼ), norm(Rₖ)
    ωn = Rᵢ ⋅ cross(Rⱼ, Rₖ)
    ωd = rᵢ*rⱼ*rₖ + rᵢ*(Rⱼ⋅Rₖ) + rⱼ*(Rₖ⋅Rᵢ) + rₖ*(Rᵢ⋅Rⱼ)
    ω = 2*atan(ωn/ωd)

    return Rᵢ, F, ω
end

"""
    compute_potential(p::PolyhedronGravity{T}, pos::AbstractVector{<:Number}, 
        G, ρ, args...) where T

Compute the gravitational potential at a given position due to a polyhedron.

# Arguments
- `p`: `PolyhedronGravity` object representing the polyhedron.
- `pos`: Position vector where the potential is computed.
- `G`: Gravitational constant.
- `ρ`: Density.

"""
function compute_potential(p::PolyhedronGravity{T}, pos::AbstractVector{<:Number}, 
    G, ρ, args...; parallel=false) where T
    if parallel
        u = _compute_potential_parallel(typeof(pos[1]), p, pos)
    else
        u = _compute_potential_serial(typeof(pos[1]), p, pos)
    end
    U = 0.5 * G * ρ * u 
    return U
end

function _edge_potential_update(p, e, r)
    Rₑ, Eₑ, Lₑ = precompute_edges(p, e, r)
    return Rₑ' * Eₑ * Rₑ * Lₑ
end

function _face_potential_update(p, f, r)
    Rᵢ, F, ω = precompute_faces(p, f, r)
    return Rᵢ' * F * Rᵢ * ω
end

function _compute_potential_serial(::T, p, pos) where T
    r = SVector{3}(pos[1], pos[2], pos[3])
    ue = mapreduce(x->_edge_potential_update(p, x, r), +, p.edges)
    uf = mapreduce(x->_face_potential_update(p, x, r), +, p.faces)
    return ue - uf
end

function _compute_potential_parallel(::T, p, pos) where T
    r = SVector{3}(pos[1], pos[2], pos[3])
    
    nblocks = Threads.nthreads() 
    u_tot = Threads.Atomic{T}(0.0)
    eblock = cld(length(p.edges), nblocks)

    Threads.@threads for t in 1:nblocks 
        u_t = mapreduce(
            x->_edge_potential_update(p, x, r), +, 
            @views(p.edges[(t-1)*eblock+1:min(t*eblock, length(p.edges))])
        )
        Threads.atomic_add!(u_tot, u_t)
    end

    fblock = cld(length(p.faces), nblocks)
    Threads.@threads for t in 1:nblocks 
        u_t = mapreduce(
            x->_face_potential_update(p, x, r), +, 
            @views(p.faces[(t-1)*fblock+1:min(t*fblock, length(p.faces))])
        )
        Threads.atomic_sub!(u_tot, u_t)
    end

    return u_tot[]
end


"""
    compute_acceleration(p::PolyhedronGravity{T}, pos::AbstractVector{<:Number}, 
        G, ρ, args...) where T

Compute the gravitational compute_acceleration at a given position due to a polyhedron.

# Arguments
- `p`: `PolyhedronGravity` object representing the polyhedron.
- `pos`: Position vector where the potential is computed.
- `G`: Gravitational constant.
- `ρ`: Density.
"""
function compute_acceleration(p::PolyhedronGravity{T}, pos::AbstractVector{<:Number}, 
    G, ρ, args...) where T
    r = SVector{3}(pos[1], pos[2], pos[3])

    δue = SVector{3, T}(0, 0, 0)
    # Loop over edges 
    for e in p.edges    
        Rₑ, Eₑ, Lₑ = precompute_edges(p, e, r)
        δue += Eₑ * Rₑ * Lₑ
    end

    δuf = SVector{3, T}(0, 0, 0)
    # Loop over faces
    for f in p.faces 
        Rᵢ, F, ω = precompute_faces(p, f, r)
        δuf += F * Rᵢ * ω
    end

    U = G * ρ * (- δue + δuf)
    return U
end

