export GravityPolyhedronData

"""

    GravityPolyhedronData{T} <: AbstractGravityModelData

A `GravityPolyhedronData` instance represents data for a polyhedron used in gravity modeling. 
It contains information about vertices, faces, and adjacency relationships between faces.

# Fields
- `vertices`: Vector of 3D vectors representing vertices.
- `faces`: Vector of 3-element vectors representing faces topology.
- `adj`: Dictionary storing the faces adjacent to each unique edge.

--- 

    GravityPolyhedronData{T}(filename::AbstractString)

Create a `GravityPolyhedronData` object from an `.obj` file.
This format requires vertices and faces to be in the same file and be organized using a prefix:
`v` for vertices and `f` for faces. Lines starting with `#` are comments and will be ingnored.

---

    GravityPolyhedronData{T}(nodefile::AbstractString, facefile::AbstractString)

Create a `GravityPolyhedronData` object from `.node` and `.face` files.
This format require a separated file for nodes and faces. Each file contains three columns with 
data. Lines starting with `#` are comments and will be ingnored.

---

    GravityPolyhedronData{T}(vertices, faces)

Create a `GravityPolyhedronData` object from vertices and faces.

"""
struct GravityPolyhedronData{T} <: AbstractGravityModelData
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

function _parse_row(::Type{T}, parts, first) where T 
    return SVector{3, T}(
        parse(T, parts[first]), parse(T, parts[first+1]), parse(T, parts[first+2])
    )
end


"""
    parse_objfile(::Type{T}, filename::AbstractString)

Parses an `.obj` file and extracts vertex and face data.

Returs `vertices`, that is an array of `SVector{3, T}` containing vertex coordinates,
and `faces` that is an array of `SVector{3, Int}` representing faces as vertex indices.
"""
function parse_objfile(::Type{T}, filename::AbstractString) where T
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
                    _parse_row(T, parts, 2)
                )
            elseif parts[1] == "f" && length(parts) == 4
                push!(
                    faces, 
                    _parse_row(Int, parts, 2)
                )
            else 
                throw( ErrorException("Line '$line' cannot be parsed!") )
            end
        end
    end
    return vertices, faces
end

"""
    parse_nodefile(::Type{T}, filename::AbstractString)

Parses an `.node` or `.face` files and extracts vertex or face data.
Returs `data`, that is an array of `SVector{3, T}` containing vertex or faces coordinates.
"""
function parse_nodefile(::Type{T}, filename::AbstractString) where T 
    data = SVector{3, T}[]
    open(filename, "r") do file 
        for line in eachline(file)
            parts = split(strip(line))
            if isempty(parts) || line[1] == '#'
                @debug "Comment line '$line' skipped"
                continue
            elseif length(parts) == 3
                push!(
                    data, 
                    _parse_row(T, parts, 1)
                )
            else 
                throw( ErrorException("Line '$line' cannot be parsed!") )
            end
        end
    end
    return data
end

function GravityPolyhedronData{T}(filename::AbstractString) where T
    if !endswith(filename, ".obj")
        throw(ErrorException("File '$filename' cannot be processed, unknown format."))
    end
    vertices, faces = parse_objfile(T, filename)
    return GravityPolyhedronData{T}(vertices, faces)
end

function GravityPolyhedronData{T}(nodefile::AbstractString, facefile::AbstractString) where T 
    if !endswith(nodefile, ".node")
        throw(ErrorException("File '$nodefile' cannot be processed, unknown format, shall be '.node'."))
    end
    if !endswith(facefile, ".face")
        throw(ErrorException("File '$facefile' cannot be processed, unknown format, shall be '.face'."))
    end

    vertices = parse_nodefile(T, nodefile)
    faces = parse_nodefile(Int, facefile)
    return GravityPolyhedronData{T}(vertices, faces)
end

function  GravityPolyhedronData{T}(vertices, faces) where T 
    adj = unique_edges(faces)
    return GravityPolyhedronData{T}(vertices, faces, adj)
end

function parse_data(::Type{T}, ::Type{GravityPolyhedronData}, filename::AbstractString) where T
    return GravityPolyhedronData{T}(filename)
end

@inline vertices(p::GravityPolyhedronData) = p.vertices
@inline faces(p::GravityPolyhedronData) = p.faces
@inline face(p::GravityPolyhedronData, i) = @inbounds faces(p)[i]
@inline edges(p::GravityPolyhedronData) = p.adj

# triplet of vertices associated to a face
@inbounds function triplet(p::GravityPolyhedronData{T}, i) where T
    f = face(p, i)
    return vertices(p)[f[1]], vertices(p)[f[2]], vertices(p)[f[3]]
end

# compute the normal of a face 
function face_normal(p::GravityPolyhedronData{T}, i) where T 
    p1, p2, p3 = triplet(p, i)
    e1 = p2 - p1 
    e2 = p3 - p2
    return unitcross(e1, e2)
end

# compute the edges vectors of a face
function face_edges(p::GravityPolyhedronData{T}, i) where T
    p1, p2, p3 = triplet(p, i)
    e1 = p2 - p1 
    e2 = p3 - p2 
    e3 = p1 - p3 
    return e1, e2, e3
end

# return the i face edges indexes (properly ordered)
function face_edges_index(p::GravityPolyhedronData{T}, i) where T 
    @inbounds v1, v2, v3 = faces(p)[i]
    return (v1, v2), (v2, v3), (v3, v1)
end

# find edge indexes given an unordered set of indexes
function find_edge(p::GravityPolyhedronData{T}, f, edge_idxs) where T
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
function face_centroid(p::GravityPolyhedronData{T}, i) where T
    p1, p2, p3 = triplet(p, i)
    return (p1 + p2 + p3) / 3
end

"""
    FaceProperties{T}

Stores properties of a face in a polyhedron gravity model.

### Fields
- `vertices`: Vertices of the face.
- `normal`: Normal vector of the face.
- `dyad`: Dyad formed by the normal vector.

--- 

    FaceProperties(p::GravityPolyhedronData{T}, i) where T

Constructs a `FaceProperties` object from a [`GravityPolyhedronData`](@ref) object and 
an index `i` representing the face.
"""
struct FaceProperties{T}
    vertices::SVector{3, Int}
    normal::SVector{3, T}
    dyad::SMatrix{3, 3, T, 9}
end

function FaceProperties(p::GravityPolyhedronData{T}, i) where T
    @inbounds v = faces(p)[i]
    n = face_normal(p, i)
    
    # dyads
    F = n * n'
    return FaceProperties{T}(v, n, F)
end

"""
    EdgeProperties{T}

Stores properties of an edge in a polyhedron gravity model.

### Fields
- `vertices`: Vertices of the edge.
- `magnitude`: Length of the edge.
- `dyad`: Dyad formed by the edge and the normals associated to its adjacent faces.

--- 

    EdgeProperties(p::GravityPolyhedronData{T}, adjfaces, faceprops) where T

Constructs an `EdgeProperties` object from a `GravityPolyhedronData`, `adjfaces`, 
containing the IDs of the two adjacent faces associated with the edge, 
and `faceprops`, a vector containing all the [`FaceProperties`](@ref) of the model.
"""
struct EdgeProperties{T}
    vertices::Set{Int}
    magnitude::T
    dyad::SMatrix{3, 3, T, 9}
end

function EdgeProperties(p::GravityPolyhedronData{T}, adjfaces, fprops) where T
    @inbounds f₁, f₂ = edges(p)[adjfaces]
    
    # indexes of the two edges for the respective faces
    # find the first one, the one associated to the adjacent face is simply the opposite 
    v₁, v₂ = find_edge(p, f₁, adjfaces)

    # edges
    e₁ = vertices(p)[v₂] - vertices(p)[v₁]

    # normals
    n₁, n₂ = fprops[f₁].normal, fprops[f₂].normal 
    ne₁, ne₂ = unitcross(e₁, n₁), unitcross(-e₁, n₂)

    # dyads
    E = n₁*ne₁' + n₂*ne₂'
    return EdgeProperties{T}(Set(adjfaces), norm(e₁), E)
end
