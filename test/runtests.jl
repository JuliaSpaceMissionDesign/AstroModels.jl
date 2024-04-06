using Test

@time begin 
    @testset "Gravity" verbose=true begin
        include("Gravity/Gravity.jl") 
    end
end;