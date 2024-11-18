export GravityPolyhedron

"""
    GravityPolyhedron{T} <: AbstractGravityModel

Efficiently handles polyhedron gravity models. 
It stores a reorganized topology of the polyhedron. 
In `vertices`, it stores the locations of the vertices of the model; in `faces`, it stores 
a collection of [`FaceProperties`](@ref), one for each face; and in `edges`, it stores a collection 
of [`EdgeProperties`](@ref), one for each unique edge.

### Fields
- `vertices`: Locations of the vertices of the model.
- `faces`: Collection of `FaceProperties`, one for each face.
- `edges`: Collection of `EdgeProperties`, one for each unique edge.

---

    GravityPolyhedron(p::GravityPolyhedronData{T}) where T

Constructs a `GravityPolyhedron` type from a [`GravityPolyhedronData`](@ref) object.
    
### References 
- Werner, R. A., & Scheeres, D. J. (1996). Exterior gravitation of a polyhedron derived and 
  compared with harmonic and mascon gravitation representations of asteroid 4769 Castalia. 
  Celestial Mechanics and Dynamical Astronomy, 65, 313-344.
"""
struct GravityPolyhedron{T} <: AbstractGravityModel{T}
    vertices::Vector{SVector{3,T}}
    faces::Vector{FaceProperties{T}}
    edges::Vector{EdgeProperties{T}}
end

function GravityPolyhedron{N}(p::GravityPolyhedronData{T}) where {N,T}
    faces_prop = Vector{FaceProperties{T}}()
    for f in eachindex(faces(p))
        push!(
            faces_prop,
            FaceProperties(p, f)
        )
    end
    edges_prop = Vector{EdgeProperties{T}}()
    for e in keys(edges(p))
        push!(
            edges_prop,
            EdgeProperties(p, e, faces_prop)
        )
    end
    return GravityPolyhedron{N}(p.vertices, faces_prop, edges_prop)
end

GravityPolyhedron(p::GravityPolyhedronData{T}) where {T} = GravityPolyhedron{T}(p)

function parse_model(::Type{T}, ::Type{GravityPolyhedron}, data::GravityPolyhedronData{N},
    args...) where {N,T}
    return GravityPolyhedron{T}(data)
end

function parse_model(::Type{T}, ::Type{GravityPolyhedron}, filename::AbstractString,
    args...) where {T}
    data = GravityPolyhedronData{T}(filename)
    return GravityPolyhedron{T}(data)
end