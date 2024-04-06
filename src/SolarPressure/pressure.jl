abstract type AbstractSunPressureModel <: AbstractJSMDModel end

struct InverseSquareSunPressure <: AbstractSunPressureModel end 

const G_SC = 1_367.0 # Solar constant at 1 astronomical unit (AU) in W/m^2
const SPEED_OF_LIGHT = 299_792_458.0 # Speed of light in vacuum in m/s
const P₀ = G_SC / SPEED_OF_LIGHT * 1e-3 # Solar pressure constant at 1 AU in mPa

"""
    compute_srp_magnitude(S::Number)

Compute solar pressure for a given distance `S` from the Sun. 

!!! warning "Units constraint"

    The solar pressure is returned in `mPa` e.g. to be used with `m^2` areas and `km` distances.
"""
@inline function compute_srp_pressure(::InverseSquareSunPressure, s::Number)
    tmp = 149_597_870.7 / s # Calculate the distance ratio (AU divided by distance from the Sun)
    return P₀ * tmp^2 
end

const INVERSE_SQUARE_SUN_PRESSURE = InverseSquareSunPressure()