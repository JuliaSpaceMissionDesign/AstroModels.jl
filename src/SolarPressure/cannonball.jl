"""
  CannonballSrpData{T}

A mutable struct representing data for the Cannonball solar radiation pressure model.

This struct holds the parameters needed for modeling solar radiation pressure (SRP) 
using the Cannonball model. It includes the following data:

- `Cr`: The reflectivity coefficient, typically equal to `1 + ρ`, where ρ is the
  spacecraft reflectivity.
- `Asc`: The spacecraft's equivalent area.
- `Msc`: The spacecraft's mass.

This data is used to compute SRP effects on a spacecraft using the Cannonball model.
"""
mutable struct CannonballSrpData{T} <: AbstractSolarPressureModelData{T}
  Cr::T
  Asc::T 
  Msc::T 
end

"""
    @inline function update!(d::CannonballSrpData{<:Number}, Msc::Number)

Update the mass parameter in Cannonball SRP data.
"""
@inline function update!(d::CannonballSrpData{<:Number}, Msc::Number)
  setfield!(d, :Msc, Msc)
  nothing
end

"""
  CannonballSrp{T} 

A struct representing the Cannonball solar radiation pressure model.

This struct encapsulates the data required for modeling solar radiation pressure (SRP)
using the Cannonball model. 

The Cannonball model is a simple approximation for SRP, considering constant
acceleration along the Sun-spacecraft direction, with the acceleration magnitude 
depending on the area-to-mass ratio and the reflectivity coefficient.
"""
struct CannonballSrp{T} <: AbstractSolarPressureModel{T}
  data::CannonballSrpData{T}
end

"""
    @inline update!(m::CannonballSrp{<:Number}, Msc::Number)

Update the mass parameter in a Cannonball SRP model.
"""
@inline update!(m::CannonballSrp{<:Number}, Msc::Number) = update!(m.data, Msc)

"""
    compute_srp_cannonball(ρ::T, Asc::T, Msc::T, s::AbstractVector{T}) where T 

Compute SRP acceleration with the cannonball model.

The cannonball model is the simplest way to approximate the SRP acceleration. This model 
considers the SRP acceleration to be constant along the Sun-spacecraft direction,
where the acceleration magnitude depends on the area-to-mass ratio and a 
reflectivity coefficient `Cr = 1+ρ`. `Asc` is the spacecraft equivalent are, `Msc` its 
mass, `P` the solar pressure and `s` the Sun-to-spacecraft position vector.

### References

- Zardaın, L., Farrés, A., & Puig, A. (2020). High-fidelity modeling and visualizing of 
  solar radiation pressure: a framework for high-fidelity analysis. UMBC Faculty Collection.
"""
@fastmath function compute_srp_cannonball(Cr::T, Asc::T, Msc::T, P::T, s::AbstractVector{T}) where T 
    sv = SA[s[1], s[2], s[3]]
    svn = sqrt(sv[1]*sv[1] + sv[2]*sv[2] + sv[3]*sv[3])
    svn3 = svn*svn*svn
    
    return - Cr*P*Asc/Msc * sv/svn3
end

"""
    compute_acceleration(m::CannonballSrp{T}, s::AbstractVector{T}, P::T) where T 

Compute the acceleration due to SRP using the Cannonball model with a specified solar 
pressure.

This function calculates the acceleration experienced by a spacecraft in the presence of
solar radiation pressure (SRP) using the Cannonball model. It takes a `CannonballSrp` 
model `m`, Sun-to-spacecraft position vector `s`, and the solar pressure `P`.
"""
@fastmath function compute_acceleration(m::CannonballSrp{T}, s::AbstractVector{T}, P::T) where T 
  return compute_srp_cannonball(m.data.rho, m.data.Asc, m.data.Msc, P, s)
end

"""
    compute_acceleration(m::CannonballSrp{T}, s::AbstractVector{T}, 
      ::AbstractSunPressureModel=INV_SQUARE_SRP) where T 

Compute the acceleration due to SRP using the Cannonball model with internal solar 
pressure calculation.

This function calculates the acceleration experienced by a spacecraft in the presence of
solar radiation pressure (SRP) using the Cannonball model. It takes a `CannonballSrp` 
model `m`, Sun-to-spacecraft position vector `s`, and an optional solar pressure model.
If no solar pressure model is specified, it uses the inverse-square law (`INV_SQUARE_SRP`)
for solar pressure.
"""
@fastmath function compute_acceleration(m::CannonballSrp{T}, s::AbstractVector{T}, 
  ::AbstractSunPressureModel=INV_SQUARE_SRP) where T 
  snorm = sqrt(s[1]*s[1] + s[2]*s[2] + s[3]*s[3])
  P = compute_solar_pressure(snorm)
  return compute_acceleration(m, s, P)
end