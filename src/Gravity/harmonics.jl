export GravityHarmonics

struct GravityHarmonics{N, T} <: AbstractGravityModel
    degree::Int

    # Parameters
    μ::T 
    radius::T 

    # Coefficients
    Clm::Vector{T}
    Slm::Vector{T}

    # Factors 
    η::Vector{T}

    # Cache
    cosl::Vector{T}
    sinl::Vector{T}
    Plm::Vector{T}
end

function Base.show(io::IO, ::GravityHarmonics{N, T}) where {N, T}
    println(io, "GravityHarmonics{$N, $T}()")
end

function GravityHarmonics(deg::Int, μ::T, radius::T, Clm::Vector{T}, Slm::Vector{T}) where T
    # Check if the dimension of the Clm, Slm coefficients are compatible with
    # the desired degree
    cdim = Int(1/2*(deg+1)*(deg+2))
    if length(Clm) != cdim-1 || length(Slm) != cdim-1
        throw(ErrorException("GravityHarmonics initialized with the wrong number of coefficients."))
    end

    # Pre-allocate factors 
    η = _compute_factors(T, deg)

    # Pre-initialises the cache
    cosl = zeros(T, deg+1)
    sinl = zeros(T, deg+1)
    Plm = zeros(T, Int(1/2*(deg+3)*(deg+4)))
    return GravityHarmonics{deg, T}(deg, μ, radius, Clm, Slm, η, cosl, sinl, Plm)
end

function _compute_factors(::Type{T}, deg::Int) where {T}
    η = zeros(T, Int(1/2*(deg+3)*(deg+4)))
    @inbounds η[2] = 1.
    @inbounds for j = 2:deg+2 
        k = j + 1
        crt = (k*(k+1))÷2
        @simd for k = 0:j-1 
            if k == 0
                @fastmath η[crt-j+k] = sqrt(0.5*(j+k+1)*(j-k))
            else 
                @fastmath η[crt-j+k] = sqrt((j+k+1)*(j-k))
            end
        end 
    end 
    η
end

function _parse_normalized_harmonics(d::D, degree::Int) where {D<:AbstractGravityHarmonicsData}
    if degree > d.max_degree 
        throw(DomainError(degree, "GravityHarmonics required degree ($degree) higher than the data one ($(d.max_degree))"))
    end
    cdim = Int(1/2*(degree+1)*(degree+2))
    return GravityHarmonics(degree, d.μ, d.radius, d.Clm[1:cdim-1], d.Slm[1:cdim-1])
end

function parse_model(::Type{T}, ::Type{GravityHarmonics}, d::GravityHarmonicsICGEMData{N, T}, 
    degree::Int) where {N, T}
    if !(d.norm == :fully_normalized) 
        throw(NotImplementedError("unnormalized ICGEM harmonics handling not avialable"))
    end
    return _parse_normalized_harmonics(d, degree)
end

function parse_model(::Type{T}, ::Type{GravityHarmonics}, d::GravityHarmonicsPDSData{N, T}, 
    degree::Int) where {N, T}
    if !d.normalized  
        throw(NotImplementedError("unnormalized PSD harmonics handling not avialable"))
    end
    return _parse_normalized_harmonics(d, degree)
end
