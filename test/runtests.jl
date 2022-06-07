using Test
using Random
using LinDistFlow
using HiGHS
using JuMP
using Ipopt
Random.seed!(42)


@testset "single phase 38-nodes 3 time steps" begin
    T = 3
    loadnodes = ["3", "5", "36", "9", "10", "11", "12", "13", "15", "17", "18", "19", "22", "25", 
                "27", "28", "30", "31", "32", "33", "34", "35"]

    loads = rand(length(loadnodes), T) * 1e2

    Pload = Dict(k =>     loads[indexin([k], loadnodes)[1], :] for k in loadnodes)
    Qload = Dict(k => 0.1*loads[indexin([k], loadnodes)[1], :] for k in loadnodes)

    Sbase = 1e6
    Vbase = 12.5e3

    ldf_inputs = Inputs(
        "data/singlephase38lines/master.dss", 
        "0";
        Pload=Pload, 
        Qload=Qload,
        Sbase=Sbase, 
        Vbase=Vbase, 
        v0 = 1.00,
        v_uplim = 1.05,
        v_lolim = 0.95,
        Ntimesteps = T
    );

    m = Model(HiGHS.Optimizer)
    build_ldf!(m, ldf_inputs)
    # can add objective here
    optimize!(m)

    @test termination_status(m) == MOI.OPTIMAL

end


# TODO https://github.com/lanl-ansi/PowerModelsDistribution.jl/blob/926b5acb6bf13231191500d404008cda580f240d/src/io/dss/dss_parse.jl#L807 fails in Ubuntu Actions tests

@testset "three phase openDSS" begin
    #=
    Validating against Arnold 2016 results, which requires:
    - dividing loads by 10? (according to values in Table I)
    - increasing impedances by 1.25 (Section IV, paragraph 3)
    - TODO? use equ.s (3) & (4) for voltage-dependent loads
    - allowing decisions for reactive power injection at busses 632, 675, 680, and 684
    - remove the "distributed load" on bus 670  (removing this load makes all Qvar ≈ 0 ???)
    =#
    p = Inputs(
        "data/13bus/IEEE13Nodeckt.dss", 
        "rg60";
        Pload=Dict(),
        Qload=Dict(),
        Sbase=500000, # get Sbase, Vbase from openDSS ?
        Vbase=4160, 
        v0 = 1.00,
        v_uplim = 1.05,
        v_lolim = 0.95,
        Ntimesteps = 1,
        P_up_bound=1e5,
        Q_up_bound=1e5,
        P_lo_bound=-1e5,
        Q_lo_bound=-1e5,
    );
    p.Qresource_nodes = ["632", "675", "680", "684"]
    for linecode in keys(p.Zdict)
        p.Zdict[linecode]["rmatrix"] *= 1.25
        p.Zdict[linecode]["xmatrix"] *= 1.25
    end
    delete!(p.Pload, "670")
    delete!(p.Qload, "670")
    p.Pload["633"] = p.Pload["634"]  # "ignoring the transformer between 633 and 644
    p.Qload["633"] = p.Qload["634"]  # "ignoring the transformer between 633 and 644
    delete!(p.Pload, "634")
    delete!(p.Qload, "634")

    @test typeof(p) == LinDistFlow.Inputs{LinDistFlow.ThreePhase}

    m = Model(Ipopt.Optimizer)
    build_ldf!(m, p)
    receiving_busses = collect(e[2] for e in p.edges)
    phases_into_bus = Dict(k=>v for (k,v) in zip(receiving_busses, p.phases))  # TODO add this to Inputs
    
    @objective(m, Min,
        0.5 * sum( m[:Qvar][b][phs][1]^2 for b in p.Qresource_nodes, phs in phases_into_bus[b]) +
        sum(
            (m[:vsqrd][b,phs1,1] - m[:vsqrd][b,phs2,1])^2
            for b in setdiff(p.busses, [p.substation_bus]), 
                (phs1, phs2) in [[1,2], [1,3], [2,3]] if phs1 in phases_into_bus[b] && phs2 in phases_into_bus[b]
        )
    )

    optimize!(m)
    @test termination_status(m) == MOI.LOCALLY_SOLVED

    phsa,phsb,phsc = 0,0,0
    for b in p.Qresource_nodes, phs in phases_into_bus[b], t in 1:1
        println(b, ".", phs, " ", round(value(m[:Qvar][b][phs][t])/p.Sbase, digits=4) )
        if phs == 1
            phsa += value(m[:Qvar][b][phs][t])/p.Sbase
        elseif phs == 2
            phsb += value(m[:Qvar][b][phs][t])/p.Sbase
        elseif phs == 3
            phsc += value(m[:Qvar][b][phs][t])/p.Sbase
        end
    end
    
    #=
    Not getting same values as Arnold for Qvar.
    Phase C has voltages at 0.95 in my solution but Arnold's lowest phase C voltages are near 0.97
    (Phase C is the only phase that will violate 0.95 without Qvar)
    =#

    vsqrd = get_bus_values("vsqrd", m, p)
    vmag = Dict(k => sqrt(v) for (k,v) in vsqrd)

end