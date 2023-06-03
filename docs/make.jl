using MovingBoundaryProblems1D
using Documenter

DocMeta.setdocmeta!(MovingBoundaryProblems1D, :DocTestSetup, :(using MovingBoundaryProblems1D); recursive=true)

makedocs(;
    modules=[MovingBoundaryProblems1D],
    authors="DanielVandH <danj.vandenheuvel@gmail.com> and contributors",
    repo="https://github.com/DanielVandH/MovingBoundaryProblems1D.jl/blob/{commit}{path}#{line}",
    sitename="MovingBoundaryProblems1D.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://DanielVandH.github.io/MovingBoundaryProblems1D.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/DanielVandH/MovingBoundaryProblems1D.jl",
    devbranch="main",
)
