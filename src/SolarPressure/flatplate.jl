"""
    compute_srp_flatplate(ρsi::T, ρdi::T, Ai::T, Msc::T, P::T, uni::AbstractVector{T}, 
        us::AbstractVector{T}) where T 

Compute SRP acceleration with the plat-plate model.

The N-plate model is an intermediate approximation that accounts for the variations of the SRP
acceleration due to changes on the spacecraft attitude, providing a more accurate approximation
of the acceleration. In this model the spacecraft is approximated by a collection of flat plates, 
each of them having different reflectivity properties, allowing the SRP magnitude to vary 
depending on the spacecraft's orientation with respect to the Sun-spacecraft line.

Notice that a flat plate only reflects on one side, hence if the light hits the back side of a plate its
SRP should not be counted. The only disadvantage of this approximation is that it does not 
take into account possible auto occultation between the different plates.

Here `ρsi` and `ρdi` are respectively the rates specular and diffusive reflection, `Ai` the 
plate surface, `Msc` the spacecraft mass, `P` the solar pressure at the current distance from 
the Sun, `us` the Sun-to-spacecraft unit vector and `uni` the plate normal.

### References

- Zardain, L., Farrés, A., & Puig, A. (2020). High-fidelity modeling and visualizing of 
  solar radiation pressure: a framework for high-fidelity analysis. UMBC Faculty Collection.
"""
@fastmath function compute_srp_flatplate(ρsi::T, ρdi::T, Ai::T, Msc::T, P::T,
    uni::AbstractVector{T}, us::AbstractVector{T}) where T 

    # Compute unit vectors and illumination
    cosθ = us[1]*uni[1] + us[2]*uni[2] + us[3]*uni[3]

    if cosθ < 0 
        # Compute acceleration components along sun and normal direction
        as = (1-ρsi) * us 
        an = 2*(ρsi*cosθ + ρdi/3) * uni

        # Compute total acceleration
        # Zardain, eq. 5
        return P * Ai/Msc * cosθ * ( as + an )
    else 
        return T(0.)
    end
end