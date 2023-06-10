export compute_acceleration

function _precompute_sincos!(c::AbstractVector{T}, s::AbstractVector{T}, 
    λ::T, max_deg::Int) where {T}
    @inbounds c[1], s[1] = 1., 0.
    @fastmath sl, cl = sincos(λ)
    @inbounds for i = 2:max_deg+1
        c[i] = cl*c[i-1] - sl*s[i-1]
        s[i] = cl*s[i-1] + sl*c[i-1]
    end
    nothing
end

# TODO: transfer legendre functions to JSMDUtils
function _legendre!(P::AbstractVector{T}, φ::T, deg::Int) where {T}
    @inbounds @fastmath begin 
        x, y = sincos(φ)
        P[1] = 1.0
        P[2] = sqrt(3)*x
        P[3] = sqrt(3)*y 
        for l = 2:deg+2
            k = l + 1
            crt = (k*(k+1))÷2
            prv = (k*(k-1))÷2
            pprv = ((k-1)*(k-2))÷2
            P[crt] = sqrt(1.0+0.5/l)*y*P[prv]
            for m = 0:l-1
                if m == 0
                    P[crt-l+m] = sqrt(2*l+1)/l*(sqrt(2*l-1)*x*P[prv-l+1+m] - (l-1)/sqrt(2*l-3)*P[pprv-l+2+m]);
                else 
                    P[crt-l+m] = sqrt(2*l+1)/sqrt((l+m)*(l-m))*(sqrt(2*l-1)*x*P[prv-l+1+m] - sqrt(l+m-1)*sqrt(l-m-1)/sqrt(2*l-3)*P[pprv-l+2+m]);
                end
            end

        end
    end 
    nothing 
end

# ------------------------------------------------------------------------------------------
# Potential model 
# ------------------------------------------------------------------------------------------

"""
    compute_potential(model::GravityHarmonics{N, T}, pos::AbstractVector{T}, args...) where {N, T}

Compute potential associated to a spherical harmonics gravitational model of degree `N` at 
the location specified by `pos` in the current object body-fixed frame. 
Other args are ignored.

### References 
- http://icgem.gfz-potsdam.de/str-0902-revised.pdf
"""
function compute_potential(model::GravityHarmonics{N, T}, pos::AbstractVector{T}, args...) where {N, T}
    # convert to spherical coordinates 
    r = sqrt(pos[1]*pos[1] + pos[2]*pos[2] + pos[3]*pos[3])
    λ = atan(pos[2], pos[1])
    φ = atan(pos[3]/sqrt(pos[1]^2 + pos[2]^2))

    r_ratio = model.radius/r 
    r_ratio_l = 1.

    U = 1.

    _precompute_sincos!(model.cosl, model.sinl, λ, model.degree)
    _legendre!(model.Plm, φ, model.degree)

    @inbounds for l = 0:model.degree
        U_l = 0.
        k = (l*(l+1))÷2 + 1
        
        @simd for m = 0:l
            Clm = model.Clm[k+m]
            Slm = model.Slm[k+m]
            Plm = model.Plm[k+m]

            slm = model.sinl[m+1]
            clm = model.cosl[m+1]

            U_l += Plm*(Clm*clm + Slm*slm)
        end

        U += r_ratio_l * U_l 
        r_ratio_l *= r_ratio
    end
    return U * model.μ/r 
end

# ------------------------------------------------------------------------------------------
# Acceleration model
# ------------------------------------------------------------------------------------------

"""
    compute_acceleration(model::GravityHarmonics{N, T}, pos::AbstractVector{T}, args...) where {N, T}

Compute acceleration associated to a spherical harmonics gravitational model of degree `N` at 
the location specified by `pos` in the current object body-fixed frame. 
Other args are ignored.

### References 
- http://icgem.gfz-potsdam.de/str-0902-revised.pdf
"""
function compute_acceleration(model::GravityHarmonics{N, T}, pos::AbstractVector{T}, args...) where {N, T}

    # convert to spherical coordinates 
    r = sqrt(pos[1]*pos[1] + pos[2]*pos[2] + pos[3]*pos[3])
    λ = atan(pos[2], pos[1])
    φ = atan(pos[3]/sqrt(pos[1]^2 + pos[2]^2))

    r_ratio = model.radius/r 
    r_ratio_l = r_ratio

    δUδr, δUδφ, δUδλ = 1., 0., 0.

    _precompute_sincos!(model.cosl, model.sinl, λ, model.degree)
    _legendre!(model.Plm, φ, model.degree)

    @inbounds for l = 2:model.degree 
        δUδr_l, δUδφ_l, δUδλ_l = 0., 0., 0.
        k = (l*(l+1)+1)÷2 # offset to the l degree coefficients

        @simd for m = 0:l
            Clm = model.Clm[k+m]
            Slm = model.Clm[k+m]
            Plm = model.Plm[k+m]

            slm = model.sinl[m+1]
            clm = model.cosl[m+1]

            if l == m 
                nPlm = 0. 
            else 
                nPlm = model.Plm[k+m+1]
            end 
            
            δUδr_l += (Clm*clm + Slm*slm)*Plm
            @fastmath δUδφ_l += (Clm*clm + Slm*slm)*(nPlm*model.η[k+m] - m*tan(φ)*Plm)
            δUδλ_l += m*Plm*(Slm*clm - Clm*slm)
        end
        r_ratio_l *= r_ratio

        δUδr += (l+1)*r_ratio_l*δUδr_l
        δUδφ += r_ratio_l*δUδφ_l
        δUδλ += r_ratio_l*δUδλ_l
    end

    tmp = model.μ/r 
    δUδr *= -tmp/r
    δUδφ *= tmp 
    δUδλ *= tmp

    if abs(φ) == π/2
        δa = SVector{3, T}(0.0, 0.0, (1/r)*δUδr*pos[3])
    else 
        r2 = r*r
        xy2 = pos[1]*pos[1] + pos[2]*pos[2]
        @fastmath xy = sqrt(xy2)

        δax = (1/r*δUδr - pos[3]/(r2*xy)*δUδφ)*pos[1] - 1/(xy2)*δUδλ*pos[2]
        δay = (1/r*δUδr - pos[3]/(r2*xy)*δUδφ)*pos[2] + 1/(xy2)*δUδλ*pos[1]
        δaz = 1/r*δUδr*pos[3] + xy/(r2)*δUδφ
        δa = SVector{3, T}(δax, δay, δaz)
    end
    return δa

end