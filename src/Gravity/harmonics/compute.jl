export compute_acceleration, compute_potential

"""
    precompute!(mod::GravityHarmonics{T}, pos::AbstractVector{T}, tid::Int) where T

Precompute gravitational harmonics coefficients for gravity field calculations.

This function precomputes and updates the gravitational harmonics coefficients needed for
gravity field calculations based on the spacecraft's position `pos`. It computes
intermediate values used in subsequent gravity potential and acceleration calculations.
"""
@fastmath function precompute!(mod::GravityHarmonics{T}, pos, tid::Int) where T

    # Compute sub-expressions
    x, y, z = pos[1], pos[2], pos[3]
    r = sqrt(x*x + y*y + z*z)/mod.radius
    r2inv = 1/(r*r)
    Rdr2inv = r2inv/mod.radius
    X = x*Rdr2inv
    Y = y*Rdr2inv
    Z = z*Rdr2inv

    get_tmp(mod.Vlm, pos)[1, 1, tid] = 1/r
    get_tmp(mod.Wlm, pos)[1, 1, tid] = T(0.)

    # Sectorial 
    # Montenbruck page 66, eq 3.29
    @inbounds for n = 1:mod.degree
        tmp = mod.η0[n]
        get_tmp(mod.Vlm, pos)[n+1, n+1, tid] = tmp * ( X * get_tmp(mod.Vlm, pos)[n, n, tid] - Y * get_tmp(mod.Wlm, pos)[n, n, tid] ) 
        get_tmp(mod.Wlm, pos)[n+1, n+1, tid] = tmp * ( X * get_tmp(mod.Wlm, pos)[n, n, tid] + Y * get_tmp(mod.Vlm, pos)[n, n, tid] ) 
    end

    # Zonal & tesseral 
    # Montenbruck pp 67, eq 3.30
    @inbounds for m = 0:mod.order 
        @simd for n = (m+1):(mod.degree)
            tmp1 = mod.η1[n+1, m+1] * Z 
            tmp2 = mod.η2[n+1, m+1] * r2inv 
            get_tmp(mod.Vlm, pos)[n+1, m+1, tid] = tmp1 * get_tmp(mod.Vlm, pos)[n, m+1, tid] + (n-1 >= 1 ? -tmp2 * get_tmp(mod.Vlm, pos)[n-1, m+1, tid] : 0) 
            get_tmp(mod.Wlm, pos)[n+1, m+1, tid] = tmp1 * get_tmp(mod.Wlm, pos)[n, m+1, tid] + (n-1 >= 1 ? -tmp2 * get_tmp(mod.Wlm, pos)[n-1, m+1, tid] : 0)  
        end
    end
    nothing
end

"""
    compute_potential(model::GravityHarmonics{T}, pos::AbstractVector{T}) where T

Compute the gravitational potential based on the gravity harmonics model.

This function calculates the gravitational potential at the spacecraft's position `pos`
using the gravity harmonics model specified by `model`. It relies on
precomputed coefficients using [`precompute!`](@ref).
"""
function compute_potential(mod::GravityHarmonics{T}, pos) where T
    
    tid = Threads.threadid()
    precompute!(mod, pos, tid)
    tmp = 0.0 # Initialize the result
    # See: Montenbruck page 66, eq 3.28
    # Sum from highest degree to smallest degree (minimizing numerical errors)
    for n = mod.degree:-1:0
        for m = n:-1:1
            tmp +=  get_tmp(mod.Vlm, pos)[n+1, m+1, tid] * mod.Clm[n+1, m+1] + get_tmp(mod.Wlm, pos)[n+1, m+1, tid] * mod.Slm[n+1, m+1]
        end
        tmp += get_tmp(mod.Vlm, pos)[n+1, 1] * mod.Clm[n+1, 1] # Zonal (m == 0)
    end
    return mod.μ/mod.radius * tmp
end

"""
    compute_acceleration(mod::GravityHarmonics{T}, pos::AbstractVector{T}) where T

Compute the gravitational acceleration based on the gravity harmonics model.

This function calculates the gravitational acceleration experienced by a spacecraft at
position `pos` using the gravity harmonics model specified by `mod`. It relies on
precomputed coefficients using [`precompute!`](@ref).
"""
function compute_acceleration(mod::GravityHarmonics{T}, pos::AbstractVector{T}, args...) where T 

    tid = Threads.threadid()
    precompute!(mod, pos, tid)

    ẍ = T(0.)
    ÿ = T(0.) 
    z̈ = T(0.) 

    g = mod.μ/mod.radius^2

    # Zonal 
    @inbounds @simd for n = 0:mod.degree # m = 0
        n̄ = n+1
        ẍ += mod.Clm[n̄, 1] * get_tmp(mod.Vlm, pos)[n̄+1, 2, tid] * mod.η0g[n̄, 1] # compute_η0_grad(n, 0)
        ÿ += mod.Clm[n̄, 1] * get_tmp(mod.Wlm, pos)[n̄+1, 2, tid] * mod.η0g[n̄, 1] # compute_η0_grad(n, 0)
        z̈ += mod.Clm[n̄, 1] * get_tmp(mod.Vlm, pos)[n̄+1, 1, tid] * mod.η2g[n̄, 1] #compute_η2_grad(n, 0)
    end
    
    # Sectorial and tesseral 
    @inbounds for n = 1:mod.degree 
        n̄ = n+1
        @simd for m = 1:n
            m̄ = m+1
            ẍ += ( -mod.Clm[n̄, m̄]*get_tmp(mod.Vlm, pos)[n̄+1, m̄+1, tid] -mod.Slm[n̄, m̄]*get_tmp(mod.Wlm, pos)[n̄+1, m̄+1, tid] ) * mod.η0g[n̄, m̄] +
                 ( mod.Clm[n̄, m̄]*get_tmp(mod.Vlm, pos)[n̄+1, m̄-1, tid]  +mod.Slm[n̄, m̄]*get_tmp(mod.Wlm, pos)[n̄+1, m̄-1, tid] ) * mod.η1g[n̄, m̄]
            ÿ += ( -mod.Clm[n̄, m̄]*get_tmp(mod.Wlm, pos)[n̄+1, m̄-1, tid] +mod.Slm[n̄, m̄]*get_tmp(mod.Vlm, pos)[n̄+1, m̄+1, tid] ) * mod.η0g[n̄, m̄] +
                 ( -mod.Clm[n̄, m̄]*get_tmp(mod.Wlm, pos)[n̄+1, m̄-1, tid] +mod.Slm[n̄, m̄]*get_tmp(mod.Vlm, pos)[n̄+1, m̄+1, tid] ) * mod.η1g[n̄, m̄]
            z̈ += ( mod.Clm[n̄, m̄]*get_tmp(mod.Vlm, pos)[n̄+1, m̄, tid]    +mod.Slm[n̄, m̄]*get_tmp(mod.Wlm, pos)[n̄+1, m̄, tid] ) * mod.η2g[n̄, m̄]
        end
    end

    return -g * SVector{3, T}( ẍ, ÿ, z̈ )
end