_test_factors(n, m) =  sqrt((2n+1)*factorial(n-m)/factorial(n+m)) * (m == 0 ? 1 : sqrt(2))  

function _precompute_terms!(V, W, pos, R, order)
    px, py, pz = pos[1], pos[2], pos[3]
    r² = px*px + py*py + pz*pz

    X = px * R / r²
    Y = py * R / r²
    Z = pz * R / r²
    R̄ = R*R /  r²

    # Initialize
    V[1, 1] = (R̄)^(1/2)

    # Zonal 
    for n in 1:order+1
        ηₙ₀ = (n-1)/n 
        tmp = n == 1 ? 0 : - ηₙ₀ * R̄ * V[n-1, 1]
        V[n+1, 1] = ((1 + ηₙ₀) * Z * V[n, 1] + tmp) 
    end

    # Sectorial/tesseral 
    for n in 1:order+1 
        V[n+1, n+1] = (2n-1) * (X * V[n, n] - Y * W[n, n])
        W[n+1, n+1] = (2n-1) * (X * W[n, n] + Y * V[n, n]) 
        for m in 2:n-1
            V[n+1, m] = (2n-1)/(n-m) * Z * V[n, m]
            W[n+1, m] = (2n-1)/(n-m) * Z * W[n, m]
            if n ≥ 2 
                V[n+1, m] -= (n+m-1)/(n-m) * R̄ * V[n-1, m] 
                W[n+1, m] -= (n+m-1)/(n-m) * R̄ * W[n-1, m]
            end
        end
    end

    for n in 1:order+1
        V[n+1, 1] *= sqrt(2n+1)
        Nnn = sqrt( 2(2n+1)/factorial(2n) )
        V[n+1, n+1] *= Nnn
        W[n+1, n+1] *= Nnn
        for m in 1:n-1
            Nnm = sqrt( factorial(n-m) * (2n+1) * 2 / factorial(n+m) )
            V[n+1, m+1] *= Nnm
            W[n+1, m+1] *= Nnm
        end
    end

    nothing

end

@testset "Precompute" verbose=true begin

    @testset "normalization factors" verbose=true begin     
        deg = 8
        N = zeros(deg+2, deg+2)
        Gravity.precompute_factors!(N, deg+1)
        for n = 0:deg+1 
            for m = 0:n 
                @test _test_factors(n, m) ≈ N[n+1, m+1]     
            end
        end
    end

end



@testset "ICGEM based, static model" verbose=true begin 
   
    PATH_ICGEM = artifact"XGM2016" 
    DATA_ICGEM = parse_data(
        Float64, GravityHarmonicsICGEMData, joinpath(PATH_ICGEM, "coeff.gfc"); 
        maxdegree=51
    )

    @testset "Accelerations" verbose=true begin
    
        @testset "zonal terms" verbose=true begin
            degree = 10
            mz = parse_model(Float64, GravityHarmonics, DATA_ICGEM, degree, true)
        
            for _ in 1:1000
                pos = rand(3)
                pos /= norm(pos)
                radius = 1 
                μ = 1 
        
                fd = ForwardDiff.gradient(x->compute_potential(mz, x, μ, radius), pos)
                an = compute_acceleration(mz, pos, μ, radius)
        
                @test maximum(abs.(fd-an)) ≤ 1e-14
            end
        end    
    
        @testset "all terms" verbose=true begin
            degree = 10
            m = parse_model(Float64, GravityHarmonics, DATA_ICGEM, degree, false)
        
            for _ in 1:1000
                pos = rand(3)
                pos /= norm(pos)
                radius = 1 
                μ = 1 
        
                fd = ForwardDiff.gradient(x->compute_potential(m, x, μ, radius), pos)
                an = compute_acceleration(m, pos, μ, radius)
        
                @test maximum(abs.(fd-an)) ≤ 1e-14
            end
        end
    
    end

    



end