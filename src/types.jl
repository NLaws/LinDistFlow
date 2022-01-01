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
    P_up_bound::Float64
    Q_up_bound::Float64
    P_lo_bound::Float64
    Q_lo_bound::Float64
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
    Nlte_cons=0,
    P_up_bound=1e4,
    Q_up_bound=1e4,
    P_lo_bound=-1e4,
    Q_lo_bound=-1e4,
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

    if v_lolim < 0 @error("lower voltage limit v_lolim cannot be less than zero") end
    if v_uplim < 0 @error("upper voltage limit v_uplim cannot be less than zero") end

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
        length(busses),  # Nnodes
        P_up_bound,
        Q_up_bound,
        P_lo_bound,
        Q_lo_bound,
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
        Nlte_cons=0,
        P_up_bound,
        Q_up_bound,
        P_lo_bound,
        Q_lo_bound,
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
    Nlte_cons=0,
    P_up_bound=1e4,
    Q_up_bound=1e4,
    P_lo_bound=-1e4,
    Q_lo_bound=-1e4,
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
        Nlte_cons=Nlte_cons,
        P_up_bound=P_up_bound,
        Q_up_bound=Q_up_bound,
        P_lo_bound=P_lo_bound,
        Q_lo_bound=Q_lo_bound,
    )
end


"""
    singlephase38linesInputs(;
        Pload=Dict{String, AbstractArray{Real, 1}}(), 
        Qload=Dict{String, AbstractArray{Real, 1}}(), 
        T=24,
        loadnodes = ["3", "5", "36", "9", "10", "11", "12", "13", "15", "17", "18", "19", "22", "25", 
                    "27", "28", "30", "31", "32", "33", "34", "35"],
        Sbase = 1e6,
        Vbase = 12.5e3,
        v0=1.0,
        v_uplim = 1.05,
        v_lolim = 0.95,
    )

Convenience function for creating a single phase network with 38 lines and nodes. 
Taken from:
Andrianesis et al. 2019 "Locational Marginal Value of Distributed Energy Resources as Non-Wires Alternatives"

NOTE that Inputs is a mutable struct (s.t. loads can be added later).
"""
function singlephase38linesInputs(;
    Pload=Dict{String, AbstractArray{Real, 1}}(), 
    Qload=Dict{String, AbstractArray{Real, 1}}(), 
    T=24,
    loadnodes = ["3", "5", "36", "9", "10", "11", "12", "13", "15", "17", "18", "19", "22", "25", 
                 "27", "28", "30", "31", "32", "33", "34", "35"],
    Sbase = 1e6,
    Vbase = 12.5e3,
    v0=1.0,
    v_uplim = 1.05,
    v_lolim = 0.95,
    )

    if isempty(Pload)  # fill in default loadnodes
        Pload = Dict(k => Real[] for k in loadnodes)
    end
    if isempty(Qload)  # fill in default loadnodes
        Qload = Dict(k => Real[] for k in loadnodes)
    end
    folderpath = joinpath(dirname(@__FILE__), "..", "test", "data")
    Inputs(
        joinpath(folderpath, "singlephase38lines.dss"), 
        "0", 
        joinpath(folderpath, "singlephase38linecodes.dss");
        Pload=Pload, 
        Qload=Qload,
        Sbase=Sbase, 
        Vbase=Vbase, 
        v0 = v0,
        v_uplim = v_uplim,
        v_lolim = v_lolim,
        Ntimesteps = T
    )
end