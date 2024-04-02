export PolyhedronGravityData

struct PolyhedronGravityData{T} <: AbstractGravityModelData
    vertices::Vector{SVector{3, T}}
    faces::Vector{SVector{3, Int}}
    adj::Dict{NTuple{2, Int}, Vector{Int}}
end

# find all the unique edges of the polyhedron
# return a dict having a (sorted) edge indexes as key and the adjacent faces as value 
function unique_edges(faces)
    unique_edges_adjfaces = Dict{Tuple{Int, Int}, Vector{Int}}()
    for (idx, face) in enumerate(faces)
        for i in 1:3
            edge = Tuple(sort([face[i], face[mod1(i + 1, 3)]]))
            if haskey(unique_edges_adjfaces, edge)
                push!(unique_edges_adjfaces[edge], idx)
            else
                unique_edges_adjfaces[edge] = Int[idx]
            end
        end
    end
    return unique_edges_adjfaces
end

function read_alias_waveform(::Type{T}, filename::AbstractString) where T
    vertices = SVector{3, T}[]
    faces = SVector{3, Int}[]
    open(filename, "r") do file 
        for line in eachline(file)
            parts = split(strip(line))
            if isempty(parts) || line[1] == '#'
                @debug "Comment line '$line' skipped"
                continue
            elseif parts[1] == "v" && length(parts) == 4
                push!(
                    vertices, 
                    SVector{3, T}(
                        parse(T, parts[2]), parse(T, parts[3]), parse(T, parts[4])
                    )
                )
            elseif parts[1] == "f" && length(parts) == 4
                push!(
                    faces, 
                    SVector{3, T}(
                        parse(Int, parts[2]), parse(Int, parts[3]), parse(Int, parts[4])
                    )
                )
            else 
                throw( ErrorException("Line '$line' cannot be parsed!") )
            end
        end
    end
    return vertices, faces
end

function PolyhedronGravityData{T}(filename::AbstractString) where T
    vertices, faces = read_alias_waveform(T, filename)
    adj = unique_edges(faces)
    return PolyhedronGravityData{T}(vertices, faces, adj)
end

function  PolyhedronGravityData{T}(vertices, faces) where T 
    adj = unique_edges(faces)
    return PolyhedronGravityData{T}(vertices, faces, adj)
end

function parse_data(::Type{T}, ::Type{PolyhedronGravityData}, filename::AbstractString) where T
    return PolyhedronGravityData{T}(filename)
end

@inline vertices(p::PolyhedronGravityData) = p.vertices
@inline faces(p::PolyhedronGravityData) = p.faces
@inline face(p::PolyhedronGravityData, i) = @inbounds faces(p)[i]
@inline edges(p::PolyhedronGravityData) = p.adj

# triplet of vertices associated to a face
@inbounds function triplet(p::PolyhedronGravityData{T}, i) where T
    f = face(p, i)
    return vertices(p)[f[1]], vertices(p)[f[2]], vertices(p)[f[3]]
end

# compute the normal of a face 
function face_normal(p::PolyhedronGravityData{T}, i) where T 
    p1, p2, p3 = triplet(p, i)
    e1 = p2 - p1 
    e2 = p3 - p2
    return unitcross(e1, e2)
end

# compute the edges vectors of a face
function face_edges(p::PolyhedronGravityData{T}, i) where T
    p1, p2, p3 = triplet(p, i)
    e1 = p2 - p1 
    e2 = p3 - p2 
    e3 = p1 - p3 
    return e1, e2, e3
end

# return the i face edges indexes (properly ordered)
function face_edges_index(p::PolyhedronGravityData{T}, i) where T 
    @inbounds v1, v2, v3 = faces(p)[i]
    return (v1, v2), (v2, v3), (v3, v1)
end

# find edge indexes given an unordered set of indexes
function find_edge(p::PolyhedronGravityData{T}, f, edge_idxs) where T
    e1, e2, e3 = face_edges_index(p, f)
    sedge = Set(edge_idxs)
    if Set(e1) == sedge
        return e1 
    elseif Set(e2) == sedge 
        return e2 
    else
        return e3 
    end
end

# compute centroid of the face
function face_centroid(p::PolyhedronGravityData{T}, i) where T
    p1, p2, p3 = triplet(p, i)
    return (p1 + p2 + p3) / 3
end

struct FaceProperties{T}
    vertices::SVector{3, Int}
    normal::SVector{3, T}
    dyad::SMatrix{3, 3, T, 9}
end

function FaceProperties(p::PolyhedronGravityData{T}, i) where T
    @inbounds v = faces(p)[i]
    n = face_normal(p, i)
    
    # dyads
    F = n * n'
    return FaceProperties{T}(v, n, F)
end

struct EdgeProperties{T}
    vertices::Set{Int}
    magnitude::T
    dyad::SMatrix{3, 3, T, 9}
end

function EdgeProperties(p::PolyhedronGravityData{T}, edge, fprops) where T
    @inbounds f₁, f₂ = edges(p)[edge]
    
    # indexes of the two edges for the respective faces
    # find the first one, the one associated to the adjacent face is simply the opposite 
    v₁, v₂ = find_edge(p, f₁, edge)

    # edges
    e₁ = vertices(p)[v₂] - vertices(p)[v₁]

    # normals
    n₁, n₂ = fprops[f₁].normal, fprops[f₂].normal 
    ne₁, ne₂ = unitcross(e₁, n₁), unitcross(-e₁, n₂)

    # dyads
    E = n₁*ne₁' + n₂*ne₂'
    return EdgeProperties{T}(Set(edge), norm(e₁), E)
end
