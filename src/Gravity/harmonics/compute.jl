@fastmath function precompute!(mod::GravityHarmonics{T}, pos::AbstractVector{T}) where T

    # Compute sub-expressions
    x, y, z = pos[1], pos[2], pos[3]
    r = sqrt(x*x + y*y + z*z)/mod.radius
    r2inv = 1/(r*r)
    Rdr2inv = r2inv/mod.radius
    X = x*Rdr2inv
    Y = y*Rdr2inv
    Z = z*Rdr2inv

    mod.Vlm[1, 1] = 1/r
    mod.Wlm[1, 1] = T(0.)

    # Sectorial 
    # Montenbruck page 66, eq 3.29
    @inbounds for n = 1:mod.degree
        tmp = mod.η0[n]
        mod.Vlm[n+1, n+1] = tmp * ( X * mod.Vlm[n, n] - Y * mod.Wlm[n, n] ) 
        mod.Wlm[n+1, n+1] = tmp * ( X * mod.Wlm[n, n] + Y * mod.Vlm[n, n] ) 
    end

    # Zonal & tesseral 
    # Montenbruck pp 67, eq 3.30
    @inbounds for m = 0:mod.order 
        @simd for n = (m+1):(mod.degree)
            tmp1 = mod.η1[n+1, m+1] * Z 
            tmp2 = mod.η2[n+1, m+1] * r2inv 
            mod.Vlm[n+1, m+1] = tmp1 * mod.Vlm[n, m+1] + (n-1 >= 1 ? -tmp2 * mod.Vlm[n-1, m+1] : 0) 
            mod.Wlm[n+1, m+1] = tmp1 * mod.Wlm[n, m+1] + (n-1 >= 1 ? -tmp2 * mod.Wlm[n-1, m+1] : 0)  
        end
    end
    nothing
end

function compute_potential(model::GravityHarmonics{T}, pos::AbstractVector{T}) where T
    precompute!(model, pos)
    tmp = 0.0 # Initialize the result
    # See: Montenbruck page 66, eq 3.28
    # Sum from highest degree to smallest degree (minimizing numerical errors)
    for n = model.degree:-1:0
        for m = n:-1:1
            tmp +=  model.Vlm[n+1, m+1] * model.Clm[n+1, m+1] + model.Wlm[n+1, m+1] * model.Slm[n+1, m+1]
        end
        tmp += model.Vlm[n+1, 1] * model.Clm[n+1, 1] # Zonal (m == 0)
    end
    return model.μ/model.radius * tmp
end

function compute_acceleration(mod::GravityHarmonics{T}, pos::AbstractVector{T}) where T 

    precompute!(mod, pos)

    ẍ = T(0.)
    ÿ = T(0.) 
    z̈ = T(0.) 

    g = mod.μ/mod.radius^2

    # Zonal 
    @inbounds @simd for n = 0:mod.degree # m = 0
        n̄ = n+1
        ẍ += mod.Clm[n̄, 1] * mod.Vlm[n̄+1, 2] * mod.η0g[n̄, 1] # compute_η0_grad(n, 0)
        ÿ += mod.Clm[n̄, 1] * mod.Wlm[n̄+1, 2] * mod.η0g[n̄, 1] # compute_η0_grad(n, 0)
        z̈ += mod.Clm[n̄, 1] * mod.Vlm[n̄+1, 1] * mod.η2g[n̄, 1] #compute_η2_grad(n, 0)
    end
    
    # Sectorial and tesseral 
    @inbounds for n = 1:mod.degree 
        n̄ = n+1
        @simd for m = 1:n
            m̄ = m+1
            ẍ += ( -mod.Clm[n̄, m̄]*mod.Vlm[n̄+1, m̄+1] -mod.Slm[n̄, m̄]*mod.Wlm[n̄+1, m̄+1] ) * mod.η0g[n̄, m̄] +
                 ( mod.Clm[n̄, m̄]*mod.Vlm[n̄+1, m̄-1]  +mod.Slm[n̄, m̄]*mod.Wlm[n̄+1, m̄-1] ) * mod.η1g[n̄, m̄]
            ÿ += ( -mod.Clm[n̄, m̄]*mod.Wlm[n̄+1, m̄-1] +mod.Slm[n̄, m̄]*mod.Vlm[n̄+1, m̄+1] ) * mod.η0g[n̄, m̄] +
                 ( -mod.Clm[n̄, m̄]*mod.Wlm[n̄+1, m̄-1] +mod.Slm[n̄, m̄]*mod.Vlm[n̄+1, m̄+1] ) * mod.η1g[n̄, m̄]
            z̈ += ( mod.Clm[n̄, m̄]*mod.Vlm[n̄+1, m̄]    +mod.Slm[n̄, m̄]*mod.Wlm[n̄+1, m̄] ) * mod.η2g[n̄, m̄]
        end
    end

    return -g * SVector{3, T}( ẍ, ÿ, z̈ )
end