export FlatplateSrp, FlatplateSrpData, compute_acceleration

"""

    FlatplateProperties{T}

A struct representing properties of a flat plate for SRP modeling.

This struct holds the properties of a flat plate used in modeling solar radiation
pressure (SRP). It includes the following data:

- `ρs`: The rate of specular reflection for the plate.
- `ρd`: The rate of diffusive reflection for the plate.
- `area`: The surface area of the plate.
- `n`: The unit vector normal to the plate's surface (one for each thread).

These properties are used to compute SRP effects on a spacecraft using the flat-plate
model.
"""
struct FlatplateProperties{T}
    ρs::T
    ρd::T 
    area::T
    n::Vector{Vector{T}}
end

"""

    FlatplateSrpData{T}

A mutable struct representing data for the flat-plate solar radiation pressure model.

This struct holds the parameters and properties required for modeling solar radiation
pressure (SRP) using the flat-plate model. It includes:

- `plates`: An array of `FlatplateProperties` objects, each representing a flat plate
    with specific properties.
- `Msc`: The spacecraft's mass.

This data is used to compute SRP effects on a spacecraft using the flat-plate SRP
model.
"""
struct FlatplateSrpData{T} <: AbstractSolarPressureModelData{T}
    plates::Vector{FlatplateProperties{T}}
    Msc::Vector{T} 
end

@inline function update!(d::FlatplateSrpData{<:Number}, Msc::Number)
    tid = Threads.threadid()
    d.Msc[tid] = Msc
    nothing
end

"""

    FlatplateSrp{T}

A struct representing the flat-plate solar radiation pressure model.

This struct encapsulates the data required for modeling solar radiation pressure (SRP)
using the flat-plate model. 

The flat-plate model is an intermediate approximation for SRP, considering variations
in SRP acceleration due to changes in spacecraft attitude. It approximates the
spacecraft as a collection of flat plates, each with different reflectivity
properties. This allows the SRP magnitude to vary depending on the spacecraft's
orientation relative to the Sun-spacecraft line.

This model provides a more accurate SRP approximation than the cannonball model.

See also: [`FlatplateSrpData`](@ref), [`FlatplateProperties`](@ref).
"""
struct FlatplateSrp{T} <: AbstractSolarPressureModelData{T}
    data::FlatplateSrpData{T}
end

"""

    @inline @inbounds getplate(m::FlatplateSrp{T}, i::Int) where T

Get the properties of a specific flat plate in the flat-plate SRP model.

This function retrieves the properties of a flat plate at index `i` from the `FlatplateSrp`
model `m`. Returns the `FlatplateProperties` object representing the specified flat plate.
"""
@inline getplate(m::FlatplateSrp{T}, i::Int) where T =  @inbounds begin m.data.plates[i] end

"""

    @inline getmass(m::FlatplateSrp{T}) where T

Get the spacecraft mass from the flat-plate SRP model.
"""
@inline getmass(m::FlatplateSrp{T}) where T = m.data.Msc[Threads.threadid()]

"""

    @inline update!(m::FlatplateSrp{T}, Msc::T) where T

Update the spacecraft mass parameter in a flat-plate SRP model.
"""
@inline update!(m::FlatplateSrp{T}, Msc::T) where T = update!(m.data, Msc)

"""

    compute_srp_flatplate(ρsi::T, ρdi::T, Ai::T, Msc::T, P::T, uni::AbstractVector{T}, 
        us::AbstractVector{T}) where T 

Compute SRP acceleration with the flat-plate model.

The N-plate model is an intermediate approximation that accounts for the variations of the SRP
acceleration due to changes on the spacecraft attitude, providing a more accurate approximation
of the acceleration than the cannonball model. See [`compute_srp_cannonball`](@ref).
In this model the spacecraft is approximated by a collection of flat plates, 
each of them having different reflectivity properties, allowing the SRP magnitude to vary 
depending on the spacecraft's orientation with respect to the Sun-spacecraft line.

Notice that a flat plate only reflects on one side, hence if the light hits the back side of a plate its
SRP should not be counted. The only disadvantage of this approximation is that it does not 
take into account possible auto occultation between the different plates.

Here `ρsi` and `ρdi` are respectively the rates specular and diffusive reflection, `Ai` the 
plate surface, `Msc` the spacecraft mass, `P` the solar pressure at the current distance from 
the Sun, `us` the **Sun-to-Spacecraft** unit vector and `uni` the plate normal.

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
        return T(0.) * us
    end
end

"""

    compute_acceleration(m::FlatplateSrp{T}, s::AbstractVector{T}) where T 

Compute SRP acceleration with the flat-plate model in the spacecraft frame.

The N-plate model is an intermediate approximation that accounts for the variations of the SRP
acceleration due to changes on the spacecraft attitude, providing a accurate approximation
of the acceleration. Here the acceleration is computed in the spacecraft frame, e.g. a body
fixed-frame associated to the spacecraft's center of mass.
Here `s` is the sun vector (Spacecraft-to-Sun) in the a the spacecraft body-fixed frame.

See [`FlatplateSrp`](@ref), [`FlatplateSrpData`](@ref), [`FlatplateProperties`](@ref).
"""
@fastmath function compute_acceleration(m::FlatplateSrp{T}, s::AbstractVector{T}, P::T) where T
    # Sun direction
    us = -unitvec(s)
    # Initialize results
    atot = SVector{3, T}(0., 0., 0.)
    # Get current thread 
    tid = Threads.threadid()
    for i in eachindex(m.data.plates) # Sum up the contribution of every plate
        plate = getplate(m, i) 
        atot += compute_srp_flatplate(plate.ρs, plate.ρd, plate.area, m.data.Msc[tid], P, plate.n[tid], us)
    end
    return atot
end

"""

    compute_acceleration(m::FlatplateSrp{T}, s::AbstractVector{T}, 
        ::AbstractSunPressureModel=INV_SQUARE_SRP) where T 

Compute SRP acceleration with the flat-plate model in the spacecraft frame with internal 
solar pressure computation.

This function calculates the acceleration experienced by a spacecraft in the presence of
solar radiation pressure (SRP) using the flat-plate model. It takes a `FlatplateSrp`
model `m`,  the sun vector (Spacecraft-to-Sun) in the a the spacecraft body-fixed frame`s`, 
and an optional solar pressure model.

If no solar pressure model is specified, it uses the inverse-square law (`INV_SQUARE_SRP`)
for solar pressure.
"""
@fastmath function compute_acceleration(m::FlatplateSrp{T}, s::AbstractVector{T}, 
    ::AbstractSunPressureModel=INV_SQUARE_SRP) where T
    # Compute sun pressure
    snorm = sqrt(s[1]*s[1] + s[2]*s[2] + s[3]*s[3])
    P = compute_solar_pressure(snorm)
    # Sun direction 
    us = -unitvec(s)
    # Initialize results
    atot = SVector{3, T}(0., 0., 0.)
    # Get current thread 
    tid = Threads.threadid()
    for i in eachindex(m.data.plates) # Sum up the contribution of every plate
        plate = getplate(m, i) 
        atot += compute_srp_flatplate(plate.ρs, plate.ρd, plate.area, m.data.Msc[tid], P, plate.n[tid], us)
    end
    return atot
end