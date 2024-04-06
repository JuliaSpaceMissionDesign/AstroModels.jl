export compute_srp_cannonball

"""
    compute_srp_cannonball(ρ, A, M, P, s::AbstractVector{<:Number}) 

Compute SRP acceleration with the cannonball model.

The cannonball model is the simplest way to approximate the SRP acceleration. This model 
considers the SRP acceleration to be constant along the Sun-spacecraft direction,
where the acceleration magnitude depends on the area-to-mass ratio and a specular reflection 
coefficient `ρ`. 
Here `A` is the spacecraft equivalent are, `M` its mass, `P` the solar pressure at the 
currenct position and `s` the Sun-to-spacecraft position vector.

### References

- Zardain, L., Farrés, A., & Puig, A. (2020). High-fidelity modeling and visualizing of 
  solar radiation pressure: a framework for high-fidelity analysis. UMBC Faculty Collection.
"""

@fastmath function _cannonball(ρ, A, M, P, sm, s)
    sm3 = sm*sm*sm
    return - (1+ρ)*P*A/M * s/sm3
end

@fastmath function compute_srp_cannonball(ρ, A, M, P, s::AbstractVector{<:Number})  
    sv = SVector{3}(s[1], s[2], s[3])
    svn = sqrt(sv[1]*sv[1] + sv[2]*sv[2] + sv[3]*sv[3])
    return _cannonball(ρ, A, M, P, svn, sv)
end

function compute_srp_cannonball(ρ, A, M, s::AbstractVector{<:Number}; pressure = INVERSE_SQUARE_SUN_PRESSURE)
    sv = SVector{3}(s[1], s[2], s[3])
    svn = sqrt(sv[1]*sv[1] + sv[2]*sv[2] + sv[3]*sv[3])
    P = compute_srp_pressure(pressure, svn)
    return _cannonball(ρ, A, M, P, svn, sv)
end