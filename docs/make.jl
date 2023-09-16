using AstrodynamicModels
using Documenter

const CI = get(ENV, "CI", "false") == "true"

makedocs(;
    modules=[AstrodynamicModels],
    authors="Andrea Pasquale <andrea.pasquale@polimi.it> and Michele Ceresoli <michele.ceresoli@polimi.it>",
    sitename="AstrodynamicModels.jl",
    format=Documenter.HTML(; prettyurls=CI, highlights=["yaml"], ansicolor=true),
    pages=[
        "Home" => "index.md", 
        "API" => 
        [
            "Gravity" => "api/Gravity.md",
            "SRP" => "api/SRP.md",
            "Atmosphere" => "api/Atmosphere.md"
        ]
    ],
    clean=true
)

# deploydocs(;
#     repo="github.com/JuliaSpaceMissionDesign/AstrodynamicModels.jl", branch="gh-pages"
# )