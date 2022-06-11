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
        joinpath("data", "singlephase38lines", "master.dss"), 
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
    - ignoring 633-634 trfx (copied line 632-633 to 633-634)
    - increasing impedances by 1.25 (Section IV, paragraph 3)
    - TODO? use equ.s (3) & (4) for voltage-dependent loads
    - allowing decisions for reactive power injection at busses 632, 675, 680, and 684
    - remove the "distributed load" on bus 670  (removing this load makes all Qvar ≈ 0 ???)
    =#
    p = Inputs(
        joinpath("data", "13bus", "IEEE13Nodeckt.dss"), 
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

    Qresource_nodes = ["632", "675", "680", "684"]

    for linecode in keys(p.Zdict)
        p.Zdict[linecode]["rmatrix"] *= 1.25
        p.Zdict[linecode]["xmatrix"] *= 1.25
    end
    delete!(p.Pload, "670")
    delete!(p.Qload, "670")
    # p.Pload["633"] = p.Pload["634"]  # "ignoring the transformer between 633 and 644
    # p.Qload["633"] = p.Qload["634"]  # "ignoring the transformer between 633 and 644
    # delete!(p.Pload, "634")
    # delete!(p.Qload, "634")

    @test typeof(p) == LinDistFlow.Inputs{LinDistFlow.ThreePhase}

    m = Model(Ipopt.Optimizer)
    build_ldf!(m, p)

    a0 = 0.85
    a1 = 0.15
    # add the var resources 
    m[:Qvar] = Dict()
    for b in Qresource_nodes
        m[:Qvar][b] = Dict()
        for phs in p.phases_into_bus[b]
            m[:Qvar][b][phs] = @variable(m, [1:p.Ntimesteps])
            delete(m, m[:cons][:injection_equalities][b][:Q][phs])
            if b in keys(p.Qload) && phs in keys(p.Qload[b])
                m[:cons][:injection_equalities][b][:Q][phs] = 
                @constraint(m, [t in 1:p.Ntimesteps],
                    m[:Qⱼ][b,phs,t] == 
                    -p.Qload[b][phs][t] / p.Sbase * (a0 + a1 * m[:vsqrd][b,phs,t]) + 
                    m[:Qvar][b][phs][t] / p.Sbase
                )
            else
                m[:cons][:injection_equalities][b][:Q][phs] = 
                @constraint(m, [t in 1:p.Ntimesteps],
                    m[:Qⱼ][b,phs,t] == 
                    m[:Qvar][b][phs][t] / p.Sbase
                )
            end
        end
    end

    
    @objective(m, Min,
        0.5* sum( (m[:Qvar][b][phs][1] / p.Sbase)^2 for b in Qresource_nodes, phs in p.phases_into_bus[b]) +
        sum(
            (m[:vsqrd][b,phs1,1] - m[:vsqrd][b,phs2,1])^2
            for b in setdiff(p.busses, [p.substation_bus]), 
                (phs1, phs2) in [[1,2], [1,3], [2,3]] if phs1 in p.phases_into_bus[b] && phs2 in p.phases_into_bus[b]
        )
    )

    # TODO looks like Arnold replaced the trfx 633-634 with a line (has 634>633 voltage on plots)
    # TODO what did Arnold replace the switch with?

    # @constraints(m, begin
    #     m[:Qvar]["632"][1][1] == 0.0081*p.Sbase*10
    #     m[:Qvar]["632"][2][1] == 0.0346*p.Sbase*10
    #     m[:Qvar]["632"][3][1] == -0.0437*p.Sbase*10
    #     m[:Qvar]["675"][1][1] == -0.005*p.Sbase*10
    #     m[:Qvar]["675"][2][1] == -0.0522*p.Sbase*10
    #     m[:Qvar]["675"][3][1] == 0.058*p.Sbase*10
    #     m[:Qvar]["680"][1][1] == -0.003*p.Sbase*10
    #     m[:Qvar]["680"][2][1] == -0.0011*p.Sbase*10
    #     m[:Qvar]["680"][3][1] == 0.0038*p.Sbase*10
    #     m[:Qvar]["684"][1][1] == 0.0003*p.Sbase*10
    #     m[:Qvar]["684"][3][1] == 0.0389*p.Sbase*10
    # end)

    optimize!(m)
    @test termination_status(m) == MOI.LOCALLY_SOLVED

    phsa,phsb,phsc = 0,0,0
    for b in Qresource_nodes, phs in p.phases_into_bus[b], t in 1:1
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
    632.1 -0.0082
    632.2 -0.0258
    632.3 0.0344
    675.1 -0.0156
    675.2 -0.0473
    675.3 0.0637
    680.1 -0.0156
    680.2 -0.0496
    680.3 0.0659
    684.1 -0.0154
    684.3 0.0631
    =#
    
    #=
    Not getting same values as Arnold for Qvar.
    Phase C has voltages at 0.95 in my solution but Arnold's lowest phase C voltages are near 0.97
    so issue is relate to min'ing the cross phase voltage differences?
    seems so: if remove control cost then get lowest voltage near 0.97
    (Phase C is the only phase that will violate 0.95 without Qvar)
    =#

    vsqrd = get_bus_values("vsqrd", m, p)
    vmag = Dict(k => sqrt(v) for (k,v) in vsqrd)

    # phs3vs = Dict()
    # for (k,v) in vmag
    #     if contains(k,"3,1]")
    #         phs3vs[k[7:9]]=v
    #     end
    # end
    # phs3vs["650"] = phs3vs["rg6"]
    # plotbusses = ["650", "632", "633", "634", "645", "646", "671", "692", "675", "680", "684", "652", "611"]
    # plotbusses3 = ["650", "632", "633", "634", "645", "646", "671", "692", "675", "680", "684", "611"] # no 652
    # plotvs = []
    # for b in plotbusses3
    #     push!(plotvs, phs3vs[b])
    # end
    # plot(plotbusses3, plotvs, line=:scatter)
    # plot!(plotbusses3, repeat([0.95], length(plotbusses3)))
    # # looks good with forced Qvar
end


@testset "modify bus injections" begin
    p = Inputs(
        joinpath("data", "13bus", "IEEE13Nodeckt.dss"), 
        "rg60";
        Pload=Dict(),
        Qload=Dict(),
        Sbase=500000,
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
    m = Model(Ipopt.Optimizer)
    build_ldf!(m, p)
    # solve with zero load on bus 680
    optimize!(m)
    phs1substationload = value.(m[:Pⱼ])[p.substation_bus, 1, 1]


    # bus 680 has three phases and no load in the IEEE13 model
    delete(m, m[:cons][:injection_equalities]["680"][:P][1])
    m[:cons][:injection_equalities]["680"][:P][1] = @constraint(m, [t in 1:p.Ntimesteps],
        m[:Pⱼ]["680",1,t] == -1e3 / p.Sbase
    )
    optimize!(m)
    @test phs1substationload < value.(m[:Pⱼ])[p.substation_bus, 1, 1]
    
end