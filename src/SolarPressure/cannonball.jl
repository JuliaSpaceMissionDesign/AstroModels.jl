"""
    compute_srp_cannonball(ρ::T, Asc::T, Msc::T, s::AbstractVector{T}) where T 

Compute SRP acceleration with the cannonball model.

The cannonball model is the simplest way to approximate the SRP acceleration. This model 
considers the SRP acceleration to be constant along the Sun-spacecraft direction,
where the acceleration magnitude depends on the area-to-mass ratio and a 
reflectivity coefficient `Cr = 1+ρ`. `Asc` is the spacecraft equivalent are, `Msc` its 
mass and `s` the Sun-to-spacecraft position vector.

### References

- Zardaın, L., Farrés, A., & Puig, A. (2020). High-fidelity modeling and visualizing of 
  solar radiation pressure: a framework for high-fidelity analysis. UMBC Faculty Collection.
"""
@fastmath function compute_srp_cannonball(ρ::T, Asc::T, Msc::T, s::AbstractVector{T}) where T 
    sv = SA[s[1], s[2], s[3]]
    svn = sqrt(sv[1]*sv[1] + sv[2]*sv[2] + sv[3]*sv[3])
    svn3 = svn*svn*svn

    # Compute solar pressure
    P = compute_solar_pressure(svn) 
    
    return - (1+ρ)*P*Asc/Msc * sv/svn3
end