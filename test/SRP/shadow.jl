# Define some test cases
test_cases = [
    (r = [1, 0, 0], s = [1, 0, 0], expected = 1),  # Satellite is directly between the body and the sun
    (r = [-1, 0, 0], s = [1, 0, 0], expected = 0),  # Satellite is directly behind the body relative to the sun
    (r = [0, 1, 0], s = [1, 0, 0], expected = 1),  # Satellite is to the side of the body relative to the sun
    (r = [0, 1, 0], s = [-1, 0, 0], expected = 1),  # Satellite is to the side of the body relative to the sun
    (r = [0, 0, 1], s = [1, 0, 0], expected = 1),  # Satellite is to the side of the body relative to the sun
]

@testset "Cylindric eclipse model" begin
    for (i, test_case) in enumerate(test_cases)
        @test SolarPressure.shadow_cylinder(test_case.r, test_case.s, 1.0) == test_case.expected
    end
end;

@testset "Conical eclipse models" begin
    Re = 6378.0
    Rs = 696340.0
    pos = [-1.5793855588225949e+03, -1.1310158524800545e+03, 6.8107001564636103e+03]
    sun = [-1.2949808021543646e+08, 6.1158901042112164e+07, 4.8494036030522384e+07]
    @test SolarPressure.shadow_conic(pos, sun, Re, Rs) ≈ 0.4888086145  atol=1e-8

    # Sunlite 
    @test SolarPressure.shadow_conic([0, 1, 0], [1e4, 0, 0], 1e-1, 1.0) ≈ 1.0
    @test SolarPressure.shadow_conic([0, 0, 1], [1e4, 0, 0], 1e-1, 1.0) ≈ 1.0
    @test SolarPressure.shadow_conic([0, -1, 0], [1e4, 0, 0], 1e-1, 1.0) ≈ 1.0
    @test SolarPressure.shadow_conic([0, 0, -1], [1e4, 0, 0], 1e-1, 1.0) ≈ 1.0

    # Occulted
    @test SolarPressure.shadow_conic([1, 0, 0], [1e4, 0, 0], 1e-1, 1.0) ≈ 0.0
end;
