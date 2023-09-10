export GravityHarmonics, parse_model

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
    Vlm::Matrix{T}
    Wlm::Matrix{T}

    η0::Vector{T}
    η1::Matrix{T}
    η2::Matrix{T}
    η0g::Matrix{T}
    η1g::Matrix{T}
    η2g::Matrix{T}
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
    Vlm = zeros(T, degree+2, order+2)
    Wlm = zeros(T, degree+2, order+2)
    η0 = zeros(T, degree+1)
    η1 = zeros(T, degree+2, order+1)
    η2 = zeros(T, degree+2, order+1)
    
    η0g = zeros(T, degree+2, order+1)
    η1g = zeros(T, degree+2, order+1)
    η2g = zeros(T, degree+2, order+1)

    for n = 1:degree+1
        η0[n] = compute_η0(n)
        η0g[n, 1] = compute_η0_grad(n-1, 0)
        η2g[n, 1] = compute_η2_grad(n-1, 0) 
    end
    for m = 0:order+1 
        for n = (m+1):degree+1
            η1[n+1, m+1] = compute_η1(n, m) 
            η2[n+1, m+1] = compute_η2(n, m) 
            η0g[n+1, m+1] = compute_η0_grad(n, m)
            η1g[n+1, m+1] = compute_η1_grad(n, m) 
            η2g[n+1, m+1] = compute_η2_grad(n, m) 
        end
    end
    return GravityHarmonics{T}(degree, order, μ, radius, Clm, Slm, Vlm, Wlm, η0, η1, η2, η0g, η1g, η2g)
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