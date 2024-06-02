using Test 

const GROUP = get(ENV, "GROUP", "All")

@time begin
    if GROUP == "All" || GROUP == "Gravity"
        include("Gravity/Gravity.jl") 
    end
    if GROUP == "All" || GROUP == "SRP"
        include("SRP/SRP.jl") 
    end
end;