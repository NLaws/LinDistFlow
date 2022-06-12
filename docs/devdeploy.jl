using LiveServer

"""
NOTE you have to dev LinDistFlow in the docs environment to get local changes
"""
function devbuildserve()
    rm("build", force=true, recursive=true)
    include("make.jl")
    cd("build")
    serve()
end

devbuildserve()