using Test

@time begin 
    @testset "Gravity" verbose=true begin
        include("Gravity/Gravity.jl") 
    end
    @testset "SRP" verbose=true begin
        include("SRP/SRP.jl") 
    end
end;