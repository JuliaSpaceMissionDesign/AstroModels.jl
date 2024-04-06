export compute_srp_flatplate

"""

    compute_srp_flatplate(ρsi, ρdi, Ai, M, P, ni::AbstractVector{<:Number}, us::AbstractVector{<:Number})

Compute SRP acceleration with the flat-plate model.

Here `ρsi` and `ρdi` are respectively the rates specular and diffusive reflection, `Ai` the 
plate surface, `M` the spacecraft mass, `P` the solar pressure at the current distance from 
the Sun, `us` the Sun-to-Spacecraft unit vector and `ni` the plate normal unit vector.

Notice that a flat plate only reflects on one side, hence if the light hits the back side of a plate its
SRP should not be counted. The only disadvantage of this approximation is that it does not 
take into account possible auto occultation between the different plates.

### References
- Montenbruck, O., Gill, E., & Lutze, F. (2002). Satellite orbits: models, methods, and  
  applications. Appl. Mech. Rev., 
- Zardain, L., Farrés, A., & Puig, A. (2020). High-fidelity modeling and visualizing of 
  solar radiation pressure: a framework for high-fidelity analysis. UMBC Faculty Collection.
"""
@fastmath function compute_srp_flatplate(ρsi, ρdi, Ai, M, P, ni::AbstractVector{<:Number}, 
    us::AbstractVector{<:Number}) 
    # Compute unit vectors and illumination
    cosθ = us[1]*ni[1] + us[2]*ni[2] + us[3]*ni[3]

    if cosθ < 0
        # Compute acceleration components along sun and normal direction
        as = (1-ρsi) * us 
        an = 2*(ρsi*cosθ + ρdi/3) * ni

        # Compute total acceleration
        # Zardain, eq. 5
        return - P * Ai/M * cosθ * ( as + an )
    else 
        return 0 * us 
    end
end

@fastmath function compute_srp_flatplate(ρsi, ρdi, Ai, M, ni::AbstractVector{<:Number}, 
    s::AbstractVector{<:Number}; pressure=INVERSE_SQUARE_SUN_PRESSURE)
    # Compute pressure
    sv = SVector{3}(s[1], s[2], s[3])
    svn = sqrt(sv[1]*sv[1] + sv[2]*sv[2] + sv[3]*sv[3])
    us = sv/svn
    P = compute_srp_pressure(pressure, svn)
    # Compute acceleration
    return compute_srp_flatplate(ρsi, ρdi, Ai, M, P, ni, us) 
end