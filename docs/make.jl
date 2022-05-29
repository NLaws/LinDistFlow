using Documenter
using LinDistFlow

makedocs(
    sitename = "LinDistFlow.jl Documentation",
    format = Documenter.HTML(),
    modules = [LinDistFlow],
    workdir = joinpath(@__DIR__, ".."),
    pages = [
        "Home" => "index.md",
        "Math" => "math.md"
    ],
)
deploydocs(
    repo = "github.com/NLaws/LinDistFlow.git",
)
