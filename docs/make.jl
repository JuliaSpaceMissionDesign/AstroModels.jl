using AstroModels
using Documenter

const CI = get(ENV, "CI", "false") == "true"

makedocs(;
    modules=[AstroModels],
    authors="Andrea Pasquale <andrea.pasquale@polimi.it>",
    sitename="AstroModels.jl",
    format=Documenter.HTML(; prettyurls=CI, highlights=["yaml"], ansicolor=true),
    pages=[
        "Home" => "index.md", 
        "API" => 
        [
            "Gravity" => "api/Gravity.md",
            # "SRP" => "api/SRP.md",
            # "Atmosphere" => "api/Atmosphere.md"
        ]
    ],
    clean=true,
)

# deploydocs(;
#     repo="github.com/JuliaSpaceMissionDesign/AstrodynamicModels.jl", branch="gh-pages"
# )