export compute_srp_cannonball


@fastmath function _cannonball(ρ, A, M, P, sm, s)
    sm3 = sm * sm * sm
    return (1 + ρ) * P * A / M * s / sm3
end

"""
    compute_srp_cannonball(ρ, A, M, P, s::AbstractVector{<:Number}) 

Compute SRP acceleration with the cannonball model.
Here `A` is the spacecraft equivalent area, `M` its mass, `P` the solar pressure at the 
currenct position and `s` the Sun-to-spacecraft position vector.

--- 

    compute_srp_cannonball(ρ, A, M, s; pressure = INVERSE_SQUARE_SUN_PRESSURE)

Compute SRP acceleration with the cannonball model, given a inverse square law for the 
Sun Pressure. 

!!! note 

    The cannonball model is the simplest way to approximate the SRP acceleration. This model 
    considers the SRP acceleration to be constant along the Sun-spacecraft direction,
    where the acceleration magnitude depends on the area-to-mass ratio and a specular reflection 
    coefficient `ρ`. 

### References
- Montenbruck, O., Gill, E., & Lutze, F. (2002). Satellite orbits: models, methods, and  
  applications. Appl. Mech. Rev., 55(2), B27-B28.
- Zardain, L., Farrés, A., & Puig, A. (2020). High-fidelity modeling and visualizing of 
  solar radiation pressure: a framework for high-fidelity analysis. UMBC Faculty Collection.
"""
@fastmath function compute_srp_cannonball(ρ, A, M, P, s::AbstractVector{<:Number})
    sv = SVector{3}(s[1], s[2], s[3])
    svn = sqrt(sv[1] * sv[1] + sv[2] * sv[2] + sv[3] * sv[3])
    return _cannonball(ρ, A, M, P, svn, sv)
end

function compute_srp_cannonball(ρ, A, M, s::AbstractVector{<:Number}; pressure=INVERSE_SQUARE_SUN_PRESSURE)
    sv = SVector{3}(s[1], s[2], s[3])
    svn = sqrt(sv[1] * sv[1] + sv[2] * sv[2] + sv[3] * sv[3])
    P = compute_srp_pressure(pressure, svn)
    return _cannonball(ρ, A, M, P, svn, sv)
end