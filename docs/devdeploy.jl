using LiveServer

function devbuildserve()
    rm("build", force=true, recursive=true)
    include("make.jl")
    cd("build")
    serve()
end

devbuildserve()