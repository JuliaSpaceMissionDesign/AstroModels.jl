export GravityHarmonics

"""
    GravityHarmonics{T} <: AbstractGravityModel

Type to efficiently handle gravity harmonics models. 
It stores both the coefficients, some precomputed coefficients as well as a computation 
cache for efficiency. Implementation according to Montenbruck (2002).

### Fields 

- `degree`: expansion degree. 
- `zonal`: boolean denoting if the computation is performed with only zonal terms. 
- `C, S`: spherical harmonics coefficients.
- `N, η0, η1, ∂η0, ∂η1, ∂η2, ∂η3`: precomputed coefficients.
- `V, W`: caches terms computation. 
    
---

    GravityHarmonics(C::AbstractMatrix{T}, S::AbstractMatrix{T}, deg::Int; 
        onlyzonal::Bool=false) 

Construct gravity harmonics cache given the coefficients, `C` and `S` and the expansion 
degree `deg`. Full expansion or zonal-only expansions could be retrieved setting `onlyzonal` 
either to `false` or `true`. 

### References
- Montenbruck, O., Gill, E., & Lutze, F. (2002). Satellite orbits: models, methods, and 
  applications. Appl. Mech. Rev., 55(2), B27-B28.
"""
struct GravityHarmonics{T} <: AbstractGravityModel{T}
    degree::Int
    zonal::Bool

    # Coefficients
    C::Matrix{T}
    S::Matrix{T}

    # Precomputed coefficients 
    N::Matrix{T}
    
    η0::Matrix{T}
    η1::Matrix{T}

    ∂η0::Vector{T}
    ∂η1::Matrix{T}
    ∂η2::Matrix{T}
    ∂η3::Matrix{T}

    # Cache
    V::DiffCache{Array{T, 3}, Vector{T}}
    W::DiffCache{Array{T, 3}, Vector{T}}
end

function Base.show(io::IO, g::GravityHarmonics{T}) where {T}
    println(io, "GravityHarmonics{$T}(degree=$(g.degree))")
end

@fastmath function precompute_factors!(N, deg)

    N[1, 1] = 1
    @inbounds for n in 1:deg 
        # m = 0
        N[n+1, 1] = sqrt(2n + 1)

        # m ≠ 0, n = m 
        Fₙ = sqrt(2n) * (2n-1) 
        N[n+1, n+1] = 1/Fₙ * N[n+1, 1] * N[n, n] * ( n==1 ? sqrt(2) : 1 )

        # m ≠ 0, n ≠ m 
        if n ≥ 2
            for m in 1:n-1
                Fₙₘ = (n+m)*(n-m+1) 
                # The following division is needed because the recurrence relation of the zonal 
                # terms have a missing factor that has to be inserted in the recurrence
                # for the general case
                m == 1 ? Fₙₘ /= 2 : nothing
                N[n+1, m+1] = 1/sqrt(Fₙₘ) * N[n+1, m]
            end
        end
    end
    nothing 
end

function precompute_coefficients!(η0, η1, deg)
    # Zonal & tesseral terms 
    η0[1, 1] = 0  
    @inbounds for n in 1:deg
        η0[n+1, 1] = (n-1)/n
        η0[n+1, n+1] = 2n-1
    end

    # Sectorial terms
    @inbounds for n in 2:deg
        for m in 1:n-1 
            η0[n+1, m+1] = (2n-1)/(n-m)
            η1[n+1, m+1] = (n+m-1)/(n-m)
        end
    end
    nothing
end

function normalization_ratio(N, n::Int, m::Int, p::Int, q::Int)
    return N[n+1, m+1]/N[p+1, q+1]
end

function factorial_ratio(n::Int, m::Int)
    # This is a simple algorithm that compute the ratio n!/m! assuming n ≥ m
    # Note that since k = n-m, this is equivalent of (m + k)!/m!, e.g. computing the 
    # product of numbers from (m + k) all the way to m
    @assert n ≥ m "factorial ratio can be computed only with n ≥ m"
    res = 1
    tmp = n
    while tmp > m
        res *= tmp
        tmp -= 1
    end
    return res
end

function normalization_ratio(p::Int, q::Int, n::Int, m::Int)
    # Computes the ratio ( CS(p, q) * VW(n, m) ) / ( CS'(p,q) * VW(n, m) )
    # where CS, VW are unnormalized coefficients and factors 
    #       CS', VW' are normalized coefficients and factors 

    @assert n ≥ 0 
    @assert m ≥ 0 
    @assert p ≥ 0 
    @assert q ≥ 0 

    num = q == 0 ? 1 : 2 
    den = m == 0 ? 1 : 2
    
    num *= 2p + 1 
    den *= 2n + 1

    if n + m > p + q 
        num *= factorial_ratio(n + m, p + q)
    else 
        den *= factorial_ratio(p + q, n + m)
    end

    if p - q > n - m 
        num *= factorial_ratio(p - q, n - m)
    else
        den *= factorial_ratio(n - m, p - q)
    end
    @fastmath return sqrt(num/den)

end

function precompute_∂coefficients!(∂η0, ∂η1, ∂η2, ∂η3, N, deg)
    @inbounds for n in 0:deg
        ∂η0[n+1] = normalization_ratio(N, n, 0, n+1, 1) 
        for m in 0:n
            if m > 0
                ∂η1[n+1, m+1] = 0.5 * normalization_ratio(N, n, m, n+1, m+1)
                ∂η2[n+1, m+1] = 0.5 * (n-m+2)*(n-m+1) * normalization_ratio(N, n, m, n+1, m-1)
            end
            ∂η3[n+1, m+1] = (n-m+1) * normalization_ratio(N, n, m, n+1, m)
        end
    end
    nothing
end

function GravityHarmonics(C::AbstractMatrix{T}, S::AbstractMatrix{T}, deg::Int; onlyzonal::Bool=false) where T 
    nth = Threads.nthreads()
    V = zeros(T, deg+2, deg+2, nth)
    W = zeros(T, deg+2, deg+2, nth)
    
    N = zeros(T, deg+2, deg+2)
    precompute_factors!(N, deg+1)

    η0 = zeros(T, deg+2, deg+2)
    η1 = zeros(T, deg+2, deg+2)
    precompute_coefficients!(η0, η1, deg+2)
    
    ∂η0 = zeros(T, deg+1)
    ∂η1 = zeros(T, deg+1, deg+1)
    ∂η2 = zeros(T, deg+1, deg+1)
    ∂η3 = zeros(T, deg+1, deg+1)
    precompute_∂coefficients!(∂η0, ∂η1, ∂η2, ∂η3, N, deg)

    return GravityHarmonics{T}(
        deg, onlyzonal, C, S, 
        N, η0, η1, ∂η0, ∂η1, ∂η2, ∂η3,
        DiffCache(V), DiffCache(W)
    )
end

"""
    terms(gh::GravityHarmonics, x)

Extracts the precomputed terms `V` and `W` from the `GravityHarmonics` object. `x` is 
necessary to access a `DiffCache` with the proper type.

This function operates in a multi-threaded environment and ensures thread safety.
"""
@inline function terms(gh::GravityHarmonics, x) 
    tid = Threads.threadid()    
    @inbounds begin
        V = @view get_tmp(gh.V, x)[:, :, tid]
        W = @view get_tmp(gh.W, x)[:, :, tid]
    end
    return V, W
end

"""
    coefficients(gh::GravityHarmonics)

Returns the coefficients `C` and `S` from the given model.
"""
@inline coefficients(gh::GravityHarmonics) = (gh.C, gh.S)

function parse_model(
    ::Type{T}, ::Type{GravityHarmonics}, data::AbstractGravityHarmonicsData{N, T},
    degree::Int, onlyzonal::Bool, args...) where {N, T}

    if !data.normalized
        throw(
            NotImplementedError("handling unnormalized coefficients currently not supported")
        )
    end

    if (typeof(data) == GravityHarmonicsICGEMData) && !(data.tide_system == :zero_tide)
        throw(
            NotImplementedError("tide system '$(d.tide_system)' not supported")
        )
    end

    if N < degree + 1 
        throw(
            ArgumentError("cannot parse GravityHarmonics of degree $(degree) with the given data")
        )
    end

    C = convert(Matrix{T}, @views(data.Clm[1:degree+1, 1:degree+1]))
    S = convert(Matrix{T}, @views(data.Slm[1:degree+1, 1:degree+1]))
    return GravityHarmonics(C, S, degree; onlyzonal=onlyzonal)
   
end
