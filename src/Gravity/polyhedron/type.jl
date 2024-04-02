export PolyhedronGravity

"""
    PolyhedronGravity{T} <: AbstractGravityModel

Type to handle the polyhedron gravity within the JSMD context.

### References 
- Werner, R. A., & Scheeres, D. J. (1996). Exterior gravitation of a polyhedron derived and 
  compared with harmonic and mascon gravitation representations of asteroid 4769 Castalia. 
  Celestial Mechanics and Dynamical Astronomy, 65, 313-344.
"""
struct PolyhedronGravity{T} <: AbstractGravityModel{T}
    vertices::Vector{SVector{3, T}}
    faces::Vector{FaceProperties{T}}
    edges::Vector{EdgeProperties{T}}
end

function PolyhedronGravity{N}(p::PolyhedronGravityData{T}) where {N, T}
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
    return PolyhedronGravity{N}(p.vertices, faces_prop, edges_prop)
end

PolyhedronGravity(p::PolyhedronGravityData{T}) where T = PolyhedronGravity{T}(p)

function parse_model(::Type{T}, ::Type{PolyhedronGravity}, data::PolyhedronGravityData{N}, 
    args...) where {N, T}
    return PolyhedronGravity{T}(data)
end

function parse_model(::Type{T}, ::Type{PolyhedronGravity}, filename::AbstractString, 
    args...) where {T}
    data = PolyhedronGravityData{T}(filename)
    return PolyhedronGravity{T}(data)
end