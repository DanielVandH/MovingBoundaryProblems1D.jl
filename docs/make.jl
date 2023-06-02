using MovingBoundaryProblem1D
using Documenter

DocMeta.setdocmeta!(MovingBoundaryProblem1D, :DocTestSetup, :(using MovingBoundaryProblem1D); recursive=true)

makedocs(;
    modules=[MovingBoundaryProblem1D],
    authors="DanielVandH <danj.vandenheuvel@gmail.com> and contributors",
    repo="https://github.com/DanielVandH/MovingBoundaryProblem1D.jl/blob/{commit}{path}#{line}",
    sitename="MovingBoundaryProblem1D.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://DanielVandH.github.io/MovingBoundaryProblem1D.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/DanielVandH/MovingBoundaryProblem1D.jl",
    devbranch="main",
)
