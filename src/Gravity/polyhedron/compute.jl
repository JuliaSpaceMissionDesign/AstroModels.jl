
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
    ωd = rᵢ * rⱼ * rₖ + rᵢ * (Rⱼ ⋅ Rₖ) + rⱼ * (Rₖ ⋅ Rᵢ) + rₖ * (Rᵢ ⋅ Rⱼ)
    ω = 2 * atan(ωn / ωd)

    return Rᵢ, F, ω
end

"""
    compute_potential(p::GravityPolyhedron{T}, pos::AbstractVector{<:Number}, 
        G, ρ, args...; parallel=false) where T

Compute the gravitational potential at a given position due to a polyhedron.

### Arguments
- `p`: `GravityPolyhedron` object representing the polyhedron.
- `pos`: Position vector where the potential is computed.
- `G`: Gravitational constant.
- `ρ`: Density.
- `parallel`: Boolean to switch on multithreaded execution. Default `false`.
"""
function compute_potential(p::GravityPolyhedron{T}, pos::AbstractVector{V},
    G, ρ, args...; parallel=false) where {T,V}
    if parallel
        u = _compute_potential_parallel(V, p, pos)
    else
        u = _compute_potential_serial(V, p, pos)
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

function _compute_potential_serial(::Type{T}, p, pos) where {T}
    r = SVector{3}(pos[1], pos[2], pos[3])
    ue = mapreduce(x -> _edge_potential_update(p, x, r), +, p.edges)
    uf = mapreduce(x -> _face_potential_update(p, x, r), +, p.faces)
    return ue - uf
end

function _compute_potential_parallel(::Type{T}, p, pos) where {T}
    r = SVector{3}(pos[1], pos[2], pos[3])
    nblocks = Threads.nthreads()
    lk = ReentrantLock()

    u_tot = zero(T)
    eblock = cld(length(p.edges), nblocks)
    Threads.@threads for t in 1:nblocks
        u_t = mapreduce(
            x -> _edge_potential_update(p, x, r), +,
            @views(p.edges[(t-1)*eblock+1:min(t * eblock, length(p.edges))])
        )
        lock(lk) do
            u_tot += u_t
        end
    end

    fblock = cld(length(p.faces), nblocks)
    Threads.@threads for t in 1:nblocks
        u_t = mapreduce(
            x -> _face_potential_update(p, x, r), +,
            @views(p.faces[(t-1)*fblock+1:min(t * fblock, length(p.faces))])
        )
        lock(lk) do
            u_tot -= u_t
        end
    end
    return u_tot
end


"""
    compute_acceleration(p::GravityPolyhedron{T}, pos::AbstractVector{<:Number}, 
        G, ρ, args...; parallel=false) where T

Compute the gravitational compute_acceleration at a given position due to a polyhedron.

### Arguments
- `p`: `GravityPolyhedron` object representing the polyhedron.
- `pos`: Position vector where the potential is computed.
- `G`: Gravitational constant.
- `ρ`: Density.
- `parallel`: Boolean to switch on multithreaded execution. Default `false`.
"""
function compute_acceleration(p::GravityPolyhedron{T}, pos::AbstractVector{V},
    G, ρ, args...; parallel=false) where {T,V}
    if parallel
        δu = _compute_acceleration_parallel(V, p, pos)
    else
        δu = _compute_acceleration_serial(V, p, pos)
    end

    U = G * ρ * δu
    return U
end

function _edge_acceleration_update(p, e, r)
    Rₑ, Eₑ, Lₑ = precompute_edges(p, e, r)
    return Eₑ * Rₑ * Lₑ
end

function _face_acceleration_update(p, f, r)
    Rᵢ, F, ω = precompute_faces(p, f, r)
    return F * Rᵢ * ω
end

function _compute_acceleration_serial(::Type{T}, p, pos) where {T}
    r = SVector{3}(pos[1], pos[2], pos[3])
    δue = mapreduce(x -> _edge_acceleration_update(p, x, r), +, p.edges)
    δuf = mapreduce(x -> _face_acceleration_update(p, x, r), +, p.faces)
    return δuf - δue
end

function _compute_acceleration_parallel(::Type{T}, p, pos) where {T}
    r = SVector{3}(pos[1], pos[2], pos[3])
    nblocks = Threads.nthreads()
    lk = ReentrantLock()

    δu_tot = SVector{3,T}(0, 0, 0)
    eblock = cld(length(p.edges), nblocks)
    Threads.@threads for t in 1:nblocks
        δu_t = mapreduce(
            x -> _edge_acceleration_update(p, x, r), +,
            @views(p.edges[(t-1)*eblock+1:min(t * eblock, length(p.edges))])
        )
        lock(lk) do
            δu_tot -= δu_t
        end
    end

    fblock = cld(length(p.faces), nblocks)
    Threads.@threads for t in 1:nblocks
        δu_t = mapreduce(
            x -> _face_acceleration_update(p, x, r), +,
            @views(p.faces[(t-1)*fblock+1:min(t * fblock, length(p.faces))])
        )
        lock(lk) do
            δu_tot += δu_t
        end
    end

    return δu_tot
end


