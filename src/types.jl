"""
    mutable struct Inputs
"""
mutable struct Inputs
    edges::Array{Tuple, 1}
    linecodes::Array{String, 1}
    linelengths::Array{Float64, 1}
    busses::Array{String}
    substation_bus::String
    Pload::Dict{String, AbstractArray{Real, 1}}
    Qload::Dict{String, AbstractArray{Real, 1}}
    Sbase::Real
    Vbase::Real
    Ibase::Real
    Zdict::Dict{String, Dict{String, Any}}
    v0::Real
    v_lolim::Real
    v_uplim::Real
    Zbase::Real
    Ntimesteps::Int
    Nequality_cons::Int
    Nlte_cons::Int
    pf::Float64
    Nnodes::Int
end


function Inputs(
    edges::Array{Tuple}, 
    linecodes::Array{String}, 
    linelengths::Array{Float64}, 
    substation_bus::String;
    Pload, 
    Qload, 
    Sbase=1, 
    Vbase=1, 
    Zdict, 
    v0, 
    v_lolim=0.95, 
    v_uplim=1.05,
    Ntimesteps=1, 
    Nequality_cons=0, 
    Nlte_cons=0
    )
    Ibase = Sbase / (Vbase * sqrt(3))
    # Ibase^2 should be used to recover amperage from lᵢⱼ ?
    Zbase = Vbase / (Ibase * sqrt(3))
    @info "Zbase: ", Zbase
    busses = String[]
    for t in edges
        push!(busses, t[1])
        push!(busses, t[2])
    end
    busses = unique(busses)

    Inputs(
        edges,
        linecodes,
        linelengths,
        busses,
        substation_bus,
        Dict(k => v/Sbase for (k,v) in Pload),
        Dict(k => v/Sbase for (k,v) in Qload),
        Sbase,
        Vbase,
        Ibase,
        Zdict,
        v0,
        v_lolim, 
        v_uplim,
        Zbase,
        Ntimesteps,
        Nequality_cons,
        Nlte_cons,
        0.1,  # power factor
        length(busses)  # Nnodes
    )
end

"""
    Inputs(
        dsslinesfilepath::String, 
        substation_bus::String, 
        dsslinecodesfilepath::String;
        Pload::AbstractDict, 
        Qload::AbstractDict, 
        Sbase=1, 
        Vbase=1, 
        v0, 
        v_lolim=0.95, 
        v_uplim=1.05,
        Ntimesteps=1, 
        Nequality_cons=0, 
        Nlte_cons=0
        )

Inputs constructor
"""
function Inputs(
    dsslinesfilepath::String, 
    substation_bus::String, 
    dsslinecodesfilepath::String;
    Pload::AbstractDict, 
    Qload::AbstractDict, 
    Sbase=1, 
    Vbase=1, 
    v0, 
    v_lolim=0.95, 
    v_uplim=1.05,
    Ntimesteps=1, 
    Nequality_cons=0, 
    Nlte_cons=0
    )
    edges, linecodes, linelengths = dss_parse_lines(dsslinesfilepath)
    linecodes_dict = dss_parse_line_codes(dsslinecodesfilepath, linecodes)
    Inputs(
        edges,
        linecodes,
        linelengths,
        substation_bus;
        Pload=Pload, 
        Qload=Qload,
        Sbase=Sbase, 
        Vbase=Vbase, 
        Zdict=linecodes_dict, 
        v0=v0,
        v_lolim = v_lolim, 
        v_uplim = v_uplim, 
        Ntimesteps=Ntimesteps,
        Nequality_cons=Nequality_cons,
        Nlte_cons=Nlte_cons
    )
end