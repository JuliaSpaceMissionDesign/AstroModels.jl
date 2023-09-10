"""
    compute_srp_flatplate(ρs::T, ρd::T, Asc::T, Msc::T, 
        s::AbstractVector{T}, n::AbstractVector{T}) where T 

Compute SRP acceleration with the plat-plate model.

The N-plate model is an intermediate approximation that accounts for the variations of the SRP
acceleration due to changes on the spacecraft attitude, providing a more accurate approximation
of the acceleration. In this model the spacecraft is approximated by a collection of flat plates, 
each of them having different reflectivity properties, allowing the SRP magnitude to vary 
depending on the spacecraft's orientation with respect to the Sun-spacecraft line.

### References

- Zardain, L., Farrés, A., & Puig, A. (2020). High-fidelity modeling and visualizing of 
  solar radiation pressure: a framework for high-fidelity analysis. UMBC Faculty Collection.
"""
@fastmath function compute_srp_flatplate(ρs::T, ρd::T, Asc::T, Msc::T, 
    s::AbstractVector{T}, n::AbstractVector{T}) where T 

    # Compute unit vectors and illumination
    us = unitvec(s)
    un = unitvec(n)
    cosθ = us[1]*un[1] + us[2]*un[2] + us[3]*un[3]

    if cosθ < 0 
        # Compute solar pressure
        P = compute_solar_pressure(svn) 
        
        # Compute acceleration components along sun and normal direction
        as = (1-ρs) * us 
        an = 2*(ρs*cosθ + ρd/3) * un

        # Compute total acceleration
        # Zardain, eq. 5
        return P * Asc/Msc * cosθ * ( as + an )
    else 
        return T(0.)
    end
end