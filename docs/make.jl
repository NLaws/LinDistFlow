using Documenter
using CommonOPF
using LinDistFlow
using JuMP

makedocs(
    sitename = "LinDistFlow.jl Documentation",
    format = Documenter.HTML(),
    modules = [LinDistFlow, CommonOPF],
    workdir = joinpath(@__DIR__, ".."),
    pages = [
        "User Documentation" => "index.md",
        "Methods" => "methods.md",
        "Math" => "math.md"
    ],
)
deploydocs(
    repo = "github.com/NLaws/LinDistFlow.git",
)
