using AstrodynamicModels
using AstrodynamicModels.Gravity
using AstrodynamicModels.SolarPressure
using AstrodynamicModels.Atmosphere

using CalcephEphemeris
using FrameTransformations

# Load ephemeris to memory
eph = load(CalcephProvider, ["/home/andrea/Documents/Kernels/spk/de440s.bsp"])
iau = load(TPC("/home/andrea/Documents/Kernels/pck/pck00011.tpc"));

# Create the graph
G = FrameSystem{2, Float64}(eph)

# Create default axes 
@axes ICRF 1 InternationalCelestialReferenceFrame
@axes IAU_EARTH 399 

# Register the ICRF axes 
add_axes_inertial!(G, ICRF)

# Register points
add_point_root!(G, :SSB, 0, ICRF)
add_point_ephemeris!(G, :EarthB, 3)
add_point_ephemeris!(G, :Sun, 10)
@point Earth 399
add_point_ephemeris!(G, Earth)
add_point_ephemeris!(G, :Moon, 301)
add_point_ephemeris!(G, :MarsB, 4)
add_point_ephemeris!(G, :JupiterB, 5)

struct NoAtmosphereModel{T} <: Atmosphere.AbstractAtmosphereModel{T} end 

struct NoSolarPressureModel{T} <: SolarPressure.AbstractSolarPressureModel{T} end

struct TwoBodyGravityModel{T} <: Gravity.AbstractGravityModel{T} end 

struct ArcModel{T,A,G,S}
    center::Symbol 
    third::Vector{Symbol}
    atmosphere::A
    gravity::G
    srp::S
end


function ArcModel{T}(center::Symbol, bodies::Vector{Symbol}, atmo::A, grav::G, srp::S) where {T, A, G, S}
    return ArcModel{T, A, G, S}(center, bodies, atmo, grav, srp)
end


BODY_MAP = Dict{Symbol, Int}(
    :Earth => 399, :Moon => 301, :Sun => 10, :SSB => 0, :EarthB => 3,
    :VenusB => 2, :MarsB => 4, :JupiterB => 5
)

GM_MAP = Dict{Symbol, Float64}(
    :SSB => NaN64, 
    :Earth => 3.9860043543609598E+05, 
    :Moon => 4.9028000661637961E+03, 
    :Sun => 1.3271244004193938E+11, 
    :EarthB => 4.0350323550225981E+05,
    :VenusB => 3.2485859200000006E+05, 
    :MarsB => 4.2828375214000022E+04,
    :JupiterB => 1.2671276480000021E+08
)

z = zeros(3,3)
sh = Gravity.GravityHarmonics{Float64}(1,2,1.,2.,z,z,z,z,zeros(3),z,z,z,z,z)

arc = ArcModel{Float64}(
    :Earth, [:Moon, :Sun, :JupiterB], 
    Atmosphere.EARTH_EXP_ATM, 
    sh, 
    NoSolarPressureModel{Float64}()
)

arc2 = ArcModel{Float64}(
    :Earth, [:Moon, :Sun, :JupiterB], 
    NoAtmosphereModel{Float64}(), 
    TwoBodyGravityModel{Float64}(), 
    NoSolarPressureModel{Float64}()
)

abstract type AbstractCelestialBody end

struct Body{T} <: AbstractCelestialBody
    id::Int
    name::Symbol
    properties::AbstractBodyProperties{T}
end

function test(arv, i)
    return arv[i].atmosphere
end