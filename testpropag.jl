using AstrodynamicModels
using AstrodynamicModels.Gravity
using AstrodynamicModels.SolarPressure
using FrameTransformations
using CalcephEphemeris

# Load ephemeris to memory
eph = load(CalcephProvider, ["/home/andrew/Documents/Kernels/spk/de440.bsp"])
iau = load(TPC("/home/andrew/Documents/Kernels/pck/pck00011.tpc"));

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

# Register body fixed frame
add_axes_bcrtod!(G, iau, Earth, IAU_EARTH, ICRF);

struct IntegratorConfig
    coi::Int 
    μc::Float64
    bodies::Vector{Int}
    μb::Vector{Float64}
    axes::Int

    frs::FrameSystem
    sh::GravityHarmonics{Float64}

    ρ::Float64
    Asc::Float64
    Msc::Float64
end

function acceleration(x, p::IntegratorConfig, j2000sec::Number)

    r = SA[x[1], x[2], x[3]]
    
    # Center gravity
    R_I2B = rotation3(G, p.axes, p.coi, j2000sec).m[1]  # Rotation matrix from ICRF to Body
    rB = R_I2B * r
    a_sh = Gravity.compute_acceleration(p.sh, rB)
    
    # Third body acceleration components
    a_3rd = SVector{3, Float64}(0., 0., 0.)
    for (i, b) in enumerate(p.bodies)
        ri = vector3(p.frs, p.coi, b, j2000sec)
        a_3rd += Gravity.compute_thirdbody_hp(p.μb[i], r, ri)
    end

    # SRP
    s = vector3(p.frs, p.coi, 10, j2000sec) + r
    sn = norm(s)
    P = SolarPressure.compute_solar_pressure(sn)
    a_srp = SolarPressure.compute_srp_cannonball(p.ρ, p.Asc, p.Bsc, P, s)

    # Total acceleration
    a_tot = R_I2B' * a_sh + a_3rd + a_srp

    return SVector{6, T}(
        x[4], x[5], x[6], a_tot[1], a_tot[2], a_tot[3]
    )

end