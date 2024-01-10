using AstroModels.Gravity

function _sph2cart(r, λ, φ)
    SA[r*cos(λ)*cos(φ), r*sin(λ)*cos(φ), r*sin(φ)]
end

function test_acceleration(model, r, λ, φ)
    pv = _sph2cart(r, λ, φ)
    return compute_acceleration(model, pv)
end

function test_potential(model, r, λ, φ)
    pv = _sph2cart(r, λ, φ)
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

@testset "gravity (w.r.t. ForwardDiff)" begin
    
    for _ in 1:10 
        λ = deg2rad(rand(0.0:0.5:360.0))
        ϕ = deg2rad(rand(-90.0:0.5:90.0))
        R = _sph2cart(XGM2016_MODEL.radius, λ, ϕ)

        ∇f_ref = ForwardDiff.gradient(x->compute_potential(XGM2016_MODEL, x), R)
        ∇f = compute_acceleration(XGM2016_MODEL, R)
        
        @test all(isapprox.(∇f_ref, ∇f; atol=1e-6, rtol=1e-8))
        @test norm((∇f_ref - ∇f)/norm(∇f)) ≤ 1e-4    
    end
end
