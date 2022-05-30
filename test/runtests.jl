using Test
using Random
using LinDistFlow
using HiGHS
using JuMP
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


@testset "three phase openDSS" begin
    p = Inputs(
        "data/13bus/IEEE13Nodeckt.dss", 
        "650";
        Pload=Dict(),  # TODO if empty load dicts try to extract from openDSS
        Qload=Dict(),
        Sbase=1, # get Sbase, Vbase from openDSS ?
        Vbase=1, 
        v0 = 1.00,
        v_uplim = 1.05,
        v_lolim = 0.95,
        Ntimesteps = 1
    )
    @test typeof(p) == LinDistFlow.Inputs{LinDistFlow.ThreePhase}
    
end