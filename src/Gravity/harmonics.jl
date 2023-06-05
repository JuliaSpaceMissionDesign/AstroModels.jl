export GravityHarmonics

"""
    GravityHarmonics{T} <: AbstractGravityModel

Type to handle the gravity harmonics models within the JSMD context. User defined `degree` 
and `order` can be supplied using the [`parse_model`](@ref) constructor associated with 
the `GravityHarmonics` type.
"""
struct GravityHarmonics{T} <: AbstractGravityModel
    degree::Int
    order::Int 

    # Parameters
    μ::T 
    radius::T 

    # Coefficients
    Clm::Matrix{T}
    Slm::Matrix{T}

    # Cache
    cosl::Vector{T}
    sinl::Vector{T}
    Plm::Matrix{T}
end

function Base.show(io::IO, g::GravityHarmonics{T}) where {T}
    println(io, "GravityHarmonics{$T}(degree=$(g.degree), order=$(g.order))")
end

function GravityHarmonics(degree::Int, order::Int, μ::T, radius::T, 
    Clm::Matrix{T}, Slm::Matrix{T}) where T

    if !all(size(Clm) .== size(Slm))
        throw(ErrorException("GravityHarmonics initialized with the wrong number of coefficients."))
    end

    if size(Clm) != (degree+1, order+1)
        throw(ErrorException("GravityHarmonics coefficients are not of dimension ($(degree+1), $(order+1))."))
    end

    cosl = zeros(T, degree+1)
    sinl = zeros(T, degree+1)
    Plm = zeros(T, degree+1, order+1)
    return GravityHarmonics{T}(degree, order, μ, radius, Clm, Slm, cosl, sinl, Plm)
end

"""
    parse_model(::Type{T}, ::Type{GravityHarmonics}, d::AbstractGravityHarmonicsData{N, T}, 
        degree::Int, order::Int=degree, args...) where {T, N}

Parse a [`GravityHarmonics`](@ref) type from a [`AbstractGravityHarmonicsData`](@ref) data container
subtype and a desired expansion `degree` and `order`.
"""
function parse_model(::Type{T}, ::Type{GravityHarmonics}, d::AbstractGravityHarmonicsData{N, T}, 
    degree::Int, order::Int=degree, args...) where {T, N}
    if d.normalized
        Clm_ = convert(Matrix{T}, @views(d.Clm[1:degree+1, 1:order+1]))
        Slm_ = convert(Matrix{T}, @views(d.Slm[1:degree+1, 1:order+1]))
        return GravityHarmonics(degree, order, d.μ, d.radius, Clm_, Slm_)
    else
        throw(NotImplementedError("handling unnormalized coefficients currently not supported"))
    end
end
