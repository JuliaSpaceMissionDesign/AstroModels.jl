export ExponentialAtmosphereData, ExponentialAtmosphere, compute_acceleration

"""
    ExponentialAtmosphereData{T}

A struct representing an exponential atmosphere model data.

This struct encapsulates the parameters defining an exponential atmosphere, which
models the variation of altitude, density, and scale height as a function of
altitude above the Earth's surface. It includes altitude levels (`h₀`), corresponding
density values (`ϱ₀`), scale heights (`H`), the Earth's radius (`radius`), and the
angular velocity (`ω`).

This struct is used to define the atmospheric properties of an exponential atmosphere
for various calculations related to atmospheric modeling and simulations.
"""
struct ExponentialAtmosphereData{T} <: AbstractAtmosphereModelData{T}
    h₀::AbstractVector{T}
    ϱ₀::AbstractVector{T}
    H::AbstractVector{T}
    radius::T 
    ω::T
end

"""
    ExponentialAtmosphere{T}

A struct representing an exponential atmosphere model.
"""
struct ExponentialAtmosphere{T} <: AbstractAtmosphereModel{T}
    data::ExponentialAtmosphereData{T}
end


const EARTH_EXP_ATM = ExponentialAtmosphere{Float64}(
    ExponentialAtmosphereData{Float64}(
        SA[0, 25, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 180, 200, 
        250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000], 
        SA[1.225, 3.899e-2, 1.774e-2, 3.972e-3, 1.057e-3, 3.206e-4, 8.770e-5, 
        1.905e-5, 3.396e-6, 5.297e-7, 9.661e-8, 2.438e-8, 8.484e-9, 3.845e-9, 
        2.070e-9, 5.464e-10, 2.789e-10, 7.248e-11, 2.418e-11, 9.518e-12, 3.725e-12, 
        1.585e-12, 6.967e-13, 1.454e-13, 3.614e-14, 1.170e-14, 5.245e-15, 3.019e-15], 
        SA[7.249, 6.349, 6.682, 7.554, 8.382, 7.714, 6.549, 5.799, 5.382, 5.877, 7.263, 
        9.473, 12.636, 16.149, 22.523, 29.740, 37.105, 45.546, 53.628, 53.298,
        58.515, 60.828, 63.822, 71.835, 88.667, 124.64, 181.05, 268.00],
        6378.1363, 
        7.292115e-5
    )
)

"""
    atmos_density_exponential(h::N, expmodel::ExponentialAtmosphere{N}) where N

Calculate the atmospheric density at a given altitude using an exponential model.
"""
function atmos_density_exponential(h::N, expmodel::ExponentialAtmosphere{N}) where N 
    # Check the bounds.
    (h < 0) && error("h must be positive.")

    # Get the values for the exponential model.
    id = (h >= 1000) ? 28 : findfirst( (x)-> x > 0 , expmodel.data.h₀ .- h ) - 1

    @inbounds begin
        h₀ = expmodel.data.h₀[id]
        ρ₀ = expmodel.data.ϱ₀[id]
        H  = expmodel.data.H[id]
    end

    # Compute the density.
    return ρ₀ * exp(-(h - h₀)/H) 
end

"""
    compute_acceleration(model::ExponentialAtmosphere{T}, pv::AbstractVector{T}, B::T) where T

Compute the acceleration experienced by an object in an exponential atmosphere.

This function calculates the acceleration experienced by an object based on its position,
velocity, and the ballistic coefficient `B` in the exponential atmosphere defined by `model`.
"""
function compute_acceleration(model::ExponentialAtmosphere{T}, 
    pv::AbstractVector{T}, B::T) where T 

    r = SA[pv[1], pv[2], pv[3]]
    v = SA[pv[4], pv[5], pv[6]]

    vrel = v - cross(SA[0, 0, model.data.ω], r)
    vr = norm(vrel)
    # it is assumed the Earth is spherical
    h = norm(r) - model.data.radius 
    ρ = atmos_density_exponential(h, model)
    # return the acceleration
    return -(0.5*ρ*B*vr) .* vrel
end

"""
    compute_acceleration(model::ExponentialAtmosphere{T}, pv::AbstractVector{T}, m::T, A::T, cD::T) where T

Compute the acceleration experienced by an object in an exponential atmosphere.

This function calculates the acceleration experienced by an object based on its position,
velocity, mass `m`, cross-sectional area `A`, and drag coefficient `cD` in the 
exponential atmosphere defined by `model`. It computes the ballistic coefficient `B` 
internally.
"""
function compute_acceleration(model::ExponentialAtmosphere{T}, 
    pv::AbstractVector{T}, m::T, A::T, cD::T) where T 
    # compute ballistic coefficient
    B = (cD*A)/m

    return compute_acceleration(model, pv, B)
end