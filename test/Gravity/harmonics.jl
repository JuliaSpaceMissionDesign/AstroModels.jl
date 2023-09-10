using AstrodynamicModels.Gravity

function test_acceleration(model, r, λ, φ)
    pv = SA[r*cos(λ)*cos(φ), r*sin(λ)*cos(φ), r*sin(φ)]
    return compute_acceleration(model, pv)
end

function test_potential(model, r, λ, φ)
    pv = SA[r*cos(λ)*cos(φ), r*sin(λ)*cos(φ), r*sin(φ)]
    return compute_potential(model, pv)
end

XGM2016_PATH = artifact"XGM2016" 
XGM2016_DATA = parse_data(Float64, GravityHarmonicsICGEMData, 
    joinpath(XGM2016_PATH, "coeff.gfc"); maxdegree=51)
XGM2016_MODEL = parse_model(Float64, GravityHarmonics, XGM2016_DATA, 50)

@testset "potential (XGM2016, spherical Earth, on surface)" begin
    raw = read_gdf_file(joinpath(XGM2016_PATH, "pot50.gdf"))

    for i in 1:size(raw)[1]
        λi = deg2rad(raw[i, 1]) 
        ϕi = deg2rad(raw[i, 2])
        value = test_potential(XGM2016_MODEL, XGM2016_MODEL.radius, λi, ϕi)
        
        @test isapprox(value, raw[i, end]/1e6, rtol=1e-14, atol=1e-12)
    end
end

@testset "gravity (XGM2016, spherical Earth, on surface)" begin
    raw = read_gdf_file(joinpath(XGM2016_PATH, "grad50.gdf"))

    for i in 1:size(raw)[1]
        λi = deg2rad(raw[i, 1]) 
        ϕi = deg2rad(raw[i, 2])
        value = norm(test_acceleration(XGM2016_MODEL, XGM2016_MODEL.radius, λi, ϕi))
        
        @test isapprox(value, raw[i, end]/1e8, rtol=1e-12, atol=1e-5)
    end
end
