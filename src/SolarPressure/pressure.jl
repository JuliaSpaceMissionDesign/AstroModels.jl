
abstract type AbstractSunPressureModel <: AbstractJSMDModel end

const SRP_AT_1AU = 1_367.0 # Solar radiation pressure at 1 astronomical unit (AU) in W/m^2
const SPEED_OF_LIGHT = 299_792_458.0 # Speed of light in vacuum in m/s
const SPR0 = SRP_AT_1AU / SPEED_OF_LIGHT # Solar pressure constant at 1 AU in N/m^2s
const AU = 149_597_870.7 # 1 astronomical unit in kilometers

struct InvSquareSunPressure <: AbstractSunPressureModel end 

const INV_SQUARE_SRP = InvSquareSunPressure()

"""
    compute_solar_pressure(S::Number)

Compute solar pressure for a given distance `S` from the Sun. Distance in kilometers.
"""
@inline function compute_solar_pressure(S::Number)
    tmp = AU / S # Calculate the distance ratio (AU divided by distance from the Sun)
    return SPR0 * tmp^2 
end