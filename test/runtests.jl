using LinDistFlow
using Test
using Random
using HiGHS
using JuMP
using Ipopt
Random.seed!(42)

# # hack for local testing
# using Pkg
# Pkg.activate("..")
# using LinDistFlow
# Pkg.activate(".")


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
    - matching the loads, capacitors, and lines in https://github.com/msankur/LinDist3Flow/tree/master/IEEE%20PES%202016/matlab/feeder13
    - changing the objective function in IEEE PES 2016/matlab/2016.01.05/V_Balance_Solver_Yapprox_20160105.m on lines 400 and 402 to Z = Z +  (FV(:,:,k1)*X)'*(FV(:,:,k1)*X); and Z = Z + 0.5*(Fu*X)'*(Fu*X);, i.e. use the squared L2 norm rather than the L2 norm 
        - note that the optimal Qvar values reported in the paper used the L2 norm objective, but the addition of the sqrt to the NLobjective does not work with Ipopt, so I modified the matlab objective from msankure to match the paper objective and the objective herein, which results in different optimal Qvar values than those reported in the paper. This test compares the optimal Qvar values to those obtained from the Matlab model with the squared L2 norm objective.
    - removing the random additions to the a0 and a1 coefficients in IEEE PES 2016/matlab/2016.01.05/runSim_iter_1step.m (and running the same file)
    - increasing impedances by 1.25 (Section IV, paragraph 3)
    - use equ.s (3) & (4) of their paper for loads (the a0 and a1 terms)
    - allowing decisions for reactive power injection at busses 632, 675, 680, and 684
    - NOTE that the signs for loads are opposite to Arnold 2016 (they assume negative injections / positive loads)
    =#

    test_values = Dict(
        "632" => [-0.003589386844896, -0.012095150336674, 0.016058346153863],
        "680" => [0.001437239423406, 0.018476869715961, -0.019964857130616],
        "675" => [0.001294902379144, 0.019607280856338, -0.021111771515487],
        "684" => [0.001114839279833, 0.0, -0.019045362995131]
    )
    p = Inputs(
        joinpath("data", "13bus", "IEEE13Nodeckt.dss"), 
        "rg60";
        Pload=Dict(),
        Qload=Dict(),
        Sbase=5_000_000, # msankur has 5,000 kvar * 1e3
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

    @test typeof(p) == LinDistFlow.Inputs{LinDistFlow.MultiPhase}

    m = Model(Ipopt.Optimizer)
    build_ldf!(m, p)

    a0 = 0.9
    a1 = 0.1

    # add loads with voltage dependency
    for b in keys(p.Pload)
        for phs in keys(p.Pload[b])
            delete(m, m[:cons][:injection_equalities][b][:P][phs])
            m[:cons][:injection_equalities][b][:P][phs] = @constraint(m, [t in 1:p.Ntimesteps],
                m[:Pj][b,phs,t] == -p.Pload[b][phs][t] / p.Sbase * (a0 + a1 * m[:vsqrd][b,phs,t])
            )
        end
    end
    for b in keys(p.Qload)
        for phs in keys(p.Qload[b])
            delete(m, m[:cons][:injection_equalities][b][:Q][phs])
            m[:cons][:injection_equalities][b][:Q][phs] = @constraint(m, [t in 1:p.Ntimesteps],
                m[:Qj][b,phs,t] == -p.Qload[b][phs][t] / p.Sbase * (a0 + a1 * m[:vsqrd][b,phs,t])
            )
        end
    end

    # add the var resources 
    m[:Qvar] = Dict()
    for b in Qresource_nodes
        m[:Qvar][b] = Dict()
        for phs in p.phases_into_bus[b]
            m[:Qvar][b][phs] = @variable(m, [1:p.Ntimesteps])
            delete(m, m[:cons][:injection_equalities][b][:Q][phs])
            if b == "675"  # capacitor injecting 200 kvar on all three phases
                m[:cons][:injection_equalities][b][:Q][phs] = 
                @constraint(m, [t in 1:p.Ntimesteps],
                    m[:Qj][b,phs,t] == 
                    -p.Qload[b][phs][t] / p.Sbase * (a0 + a1 * m[:vsqrd][b,phs,t]) + 
                    m[:Qvar][b][phs][t] / p.Sbase +
                    200_000 / p.Sbase
                )
                delete(m, m[:cons][:injection_equalities][b][:P][phs])
                m[:cons][:injection_equalities][b][:P][phs] = @constraint(m, [t in 1:p.Ntimesteps],
                    m[:Pj][b,phs,t] == -p.Pload[b][phs][t] / p.Sbase * (a0 + a1 * m[:vsqrd][b,phs,t])
                )
            elseif b in keys(p.Qload) && phs in keys(p.Qload[b])
                m[:cons][:injection_equalities][b][:Q][phs] = 
                @constraint(m, [t in 1:p.Ntimesteps],
                    m[:Qj][b,phs,t] == 
                    -p.Qload[b][phs][t] / p.Sbase * (a0 + a1 * m[:vsqrd][b,phs,t]) + 
                    m[:Qvar][b][phs][t] / p.Sbase
                )
            else
                m[:cons][:injection_equalities][b][:Q][phs] = 
                @constraint(m, [t in 1:p.Ntimesteps],
                    m[:Qj][b,phs,t] == 
                    m[:Qvar][b][phs][t] / p.Sbase
                )
            end
        end
    end

    # capacitor on bus 611, only has phase 3
    delete(m, m[:cons][:injection_equalities]["611"][:Q][3])
    m[:cons][:injection_equalities]["611"][:Q][3] = 
        @constraint(m, [t in 1:p.Ntimesteps],
            m[:Qj]["611",3,t] == 
            -p.Qload["611"][3][t] / p.Sbase * (a0 + a1 * m[:vsqrd]["611",3,t]) + 100_000 / p.Sbase
        )
    
    @objective(m, Min,
        0.5* sum( (m[:Qvar][b][phs][1] / p.Sbase)^2 for b in Qresource_nodes, phs in p.phases_into_bus[b]) +
        sum(
            (m[:vsqrd][b,phs1,1] - m[:vsqrd][b,phs2,1])^2
            for b in setdiff(p.busses, [p.substation_bus]), 
                (phs1, phs2) in [[1,2], [1,3], [2,3]] if phs1 in p.phases_into_bus[b] && phs2 in p.phases_into_bus[b]
        )
    )

    ## does not work ?!
    # @NLobjective(m, Min,
    #     0.5* sqrt(sum( (m[:Qvar][b][phs][1] / p.Sbase)^2 for b in Qresource_nodes, phs in p.phases_into_bus[b])) +
    #     sqrt(sum(
    #         (m[:vsqrd][b,phs1,1] - m[:vsqrd][b,phs2,1])^2
    #         for b in setdiff(p.busses, [p.substation_bus]), 
    #             (phs1, phs2) in [[1,2], [1,3], [2,3]] if phs1 in p.phases_into_bus[b] && phs2 in p.phases_into_bus[b]
    #     ))
    # )

    optimize!(m)
    @test termination_status(m) == MOI.LOCALLY_SOLVED
        
    for b in Qresource_nodes, phs in p.phases_into_bus[b]
        @test -test_values[b][phs] ≈ value(m[:Qvar][b][phs][1])/p.Sbase rtol=1e-3
    end
end


@testset "multiphase voltage regulator" begin

    p = Inputs(
        joinpath("data", "13bus", "IEEE13Nodeckt.dss"), 
        "rg60";
        Pload=Dict(),
        Qload=Dict(),
        Sbase=5_000_000,
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

    # only phase 3 684 -> 611
    p.regulators = Dict( ("684","611") => Dict(:turn_ratio => Dict(3 => 1.05)) )

    m = Model(HiGHS.Optimizer)
    build_ldf!(m, p)
    @objective(m, Min, sum(
        m[:Pj][p.substation_bus, phs, 1] + m[:Qj][p.substation_bus, phs, 1]
        for phs in 1:3)
    )
    optimize!(m)
    @test termination_status(m) == MOI.OPTIMAL

    vs = value.(m[:vsqrd]).data
    @test vs[("611", 3, 1)] ≈ 1.05^2 * vs[("684", 3, 1)]

    p.regulators = Dict( ("684","611") => Dict(:vreg => Dict(3 => 1.02)) )

    m = Model(HiGHS.Optimizer)
    build_ldf!(m, p)
    @objective(m, Min, sum(
        m[:Pj][p.substation_bus, phs, 1] + m[:Qj][p.substation_bus, phs, 1]
        for phs in 1:3)
    )
    optimize!(m)
    @test termination_status(m) == MOI.OPTIMAL

    vs = value.(m[:vsqrd]).data
    @test vs[("611", 3, 1)] ≈ 1.02^2

    # all three phases 632 -> 633
    p.regulators = Dict( ("632","633") => Dict(:vreg => Dict(1 => 1.03, 2 => 1.03, 3 => 1.03)) )
    m = Model(HiGHS.Optimizer)
    build_ldf!(m, p)
    @objective(m, Min, sum(
        m[:Pj][p.substation_bus, phs, 1] + m[:Qj][p.substation_bus, phs, 1]
        for phs in 1:3)
    )
    optimize!(m)
    @test termination_status(m) == MOI.OPTIMAL
    vs = value.(m[:vsqrd]).data
    @test has_vreg(p, "633")
    for phs in 1:3
        @test vs[("633", phs, 1)] ≈ 1.03^2
    end

end