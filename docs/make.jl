using MechanochemicalPatterns
using Documenter

DocMeta.setdocmeta!(MechanochemicalPatterns, :DocTestSetup, :(using MechanochemicalPatterns); recursive=true)

makedocs(;
    modules=[MechanochemicalPatterns],
    authors="Steffen Plunder <steffen.plunder@web.de> and contributors",
    repo="https://github.com/SteffenPL/MechanochemicalPatterns.jl/blob/{commit}{path}#{line}",
    sitename="MechanochemicalPatterns.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://SteffenPL.github.io/MechanochemicalPatterns.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/SteffenPL/MechanochemicalPatterns.jl",
    devbranch="main",
)
