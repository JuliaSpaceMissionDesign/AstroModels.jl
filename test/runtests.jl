using Test

@time begin 
    @time @testset "Gravity" verbose=true begin
        include("Gravity/Gravity.jl") 
    end
end;