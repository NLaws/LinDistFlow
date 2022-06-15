"""
    build_ldf!(m::JuMP.AbstractModel, p::Inputs)

Add variables and constraints to `m` using the values in `p`. Calls the following functions:
```julia
add_variables(m, p)
constrain_power_balance(m, p)
constrain_substation_voltage(m, p)
constrain_KVL(m, p)
constrain_loads(m, p)
```
"""
function build_ldf!(m::JuMP.AbstractModel, p::Inputs)

    add_variables(m, p)
    constrain_power_balance(m, p)
    constrain_substation_voltage(m, p)
    constrain_KVL(m, p)
    constrain_loads(m, p)

end


function add_variables(m, p::Inputs)
    T = 1:p.Ntimesteps
    # bus injections
    @variables m begin
        p.P_lo_bound <= Pⱼ[p.busses, T] <= p.P_up_bound
        p.Q_lo_bound <= Qⱼ[p.busses, T] <= p.Q_up_bound
    end

    # voltage squared
    @variable(m, p.v_lolim^2 <= vsqrd[p.busses, T] <= p.v_uplim^2 ) 
    
    ij_edges = [string(i*"-"*j) for j in p.busses for i in i_to_j(j, p)]

    # line flows, power sent from i to j
    @variable(m, p.P_lo_bound <= Pᵢⱼ[ij_edges, T] <= p.P_up_bound )
    
    @variable(m, p.Q_lo_bound <= Qᵢⱼ[ij_edges, T] <= p.Q_up_bound )

    nothing
end


function add_variables(m, p::Inputs{ThreePhase})
    receiving_busses = collect(e[2] for e in p.edges)
    d = Dict(k=>v for (k,v) in zip(receiving_busses, p.phases))
    d[p.substation_bus] = [1,2,3]
    T = 1:p.Ntimesteps
    # bus injections
    @variables m begin
        p.P_lo_bound <= Pⱼ[b in keys(d), d[b], T] <= p.P_up_bound
        p.Q_lo_bound <= Qⱼ[b in keys(d), d[b], T] <= p.Q_up_bound
    end
    
    # voltage squared
    @variable(m, p.v_lolim^2 <= vsqrd[b in keys(d), d[b], T] <= p.v_uplim^2 ) 
    
    ij_edges = [string(i*"-"*j) for (i,j) in p.edges]
    edge2phases = Dict(k=>v for (k,v) in zip(ij_edges, p.phases))  # this will fail for mesh network

    # line flows, power sent from i to j
    @variable(m, p.P_lo_bound <= Pᵢⱼ[e in keys(edge2phases), edge2phases[e], T] <= p.P_up_bound )
    
    @variable(m, p.Q_lo_bound <= Qᵢⱼ[e in keys(edge2phases), edge2phases[e], T] <= p.Q_up_bound )
    # TODO line flow limit inputs
    nothing
end


function constrain_power_balance(m, p::Inputs)
    Pⱼ = m[:Pⱼ]
    Qⱼ = m[:Qⱼ]
    Pᵢⱼ = m[:Pᵢⱼ]
    Qᵢⱼ = m[:Qᵢⱼ]
    # TODO change Pⱼ and Qⱼ to expressions, make P₀ and Q₀ dv's, which will reduce # of variables
    # by (Nnodes - 1)*8760 and number of constraints by 6*(Nnodes - 1)*8760
    for j in p.busses
        if isempty(i_to_j(j, p)) && !isempty(j_to_k(j, p)) # source nodes
            pcon = @constraint(m,  [t in 1:p.Ntimesteps],
                Pⱼ[j,t] - sum( Pᵢⱼ[string(j*"-"*k), t] for k in j_to_k(j, p) ) == 0
            )
            qcon = @constraint(m, [t in 1:p.Ntimesteps],
                Qⱼ[j,t] - sum( Qᵢⱼ[string(j*"-"*k), t] for k in j_to_k(j, p) ) == 0
            )
        elseif isempty(i_to_j(j, p)) && isempty(j_to_k(j, p))  # unconnected nodes
            @warn "Bus $j has no edges, setting Pⱼ and Qⱼ to zero."
            pcon = @constraint(m, [t in 1:p.Ntimesteps],
                Pⱼ[j,t] == 0
            )
            qcon = @constraint(m, [t in 1:p.Ntimesteps],
                Qⱼ[j,t] == 0
            )
        elseif !isempty(i_to_j(j, p)) && isempty(j_to_k(j, p))  # leaf nodes
            pcon = @constraint(m, [t in 1:p.Ntimesteps],
                sum( Pᵢⱼ[string(i*"-"*j), t] for i in i_to_j(j, p) ) + Pⱼ[j, t] == 0
            )
            qcon = @constraint(m, [t in 1:p.Ntimesteps],
                sum( Qᵢⱼ[string(i*"-"*j), t] for i in i_to_j(j, p) ) + Qⱼ[j, t] == 0
            )
        else
            pcon =  @constraint(m, [t in 1:p.Ntimesteps],
                sum( Pᵢⱼ[string(i*"-"*j), t] for i in i_to_j(j, p) ) +
                Pⱼ[j,t] - sum( Pᵢⱼ[string(j*"-"*k), t] for k in j_to_k(j, p) ) == 0
            )
            qcon = @constraint(m, [t in 1:p.Ntimesteps],
                sum( Qᵢⱼ[string(i*"-"*j), t] for i in i_to_j(j, p) ) +
                Qⱼ[j,t] - sum( Qᵢⱼ[string(j*"-"*k), t] for k in j_to_k(j, p) ) == 0
            )
        end
    end
    nothing
end


function constrain_power_balance(m, p::Inputs{ThreePhase})
    Pⱼ = m[:Pⱼ]
    Qⱼ = m[:Qⱼ]
    Pᵢⱼ = m[:Pᵢⱼ]
    Qᵢⱼ = m[:Qᵢⱼ]

    source_nodes = setdiff([string(e[1]) for e in p.edges],[e[2] for e in p.edges])
    for j in source_nodes
        for phs in [1,2,3]
            ks_on_phs = [k for k in j_to_k(j, p) if phs in p.phases_into_bus[k]]
            if !isempty(ks_on_phs)
                @constraint(m, [t in 1:p.Ntimesteps],
                    Pⱼ[j,phs,t] - sum( Pᵢⱼ[string(j*"-"*k), phs, t] for k in ks_on_phs ) == 0
                )
                @constraint(m, [t in 1:p.Ntimesteps],
                    Qⱼ[j,phs,t] - sum( Qᵢⱼ[string(j*"-"*k), phs, t] for k in ks_on_phs ) == 0
                )
            else
                @warn "Source node $j is not connected to any busses on phase $phs. Setting Pⱼ and Qⱼ to zero."
                @constraint(m, [t in 1:p.Ntimesteps],
                    Pⱼ[j,phs,t] == 0
                )
                @constraint(m, [t in 1:p.Ntimesteps],
                    Qⱼ[j,phs,t] == 0
                )
            end
        end
    end # source nodes

    for j in setdiff(p.busses, source_nodes)
        if isempty(i_to_j(j, p)) && isempty(j_to_k(j, p))  # unconnected nodes
            @warn "Bus $j has no edges in or out; setting Pⱼ and Qⱼ to zero."
            @constraint(m, [phs in 1:3, t in 1:p.Ntimesteps],
                Pⱼ[j,phs,t] == 0
            )
            @constraint(m, [phs in 1:3, t in 1:p.Ntimesteps],
                Qⱼ[j,phs,t] == 0
            )
        else  
            for phs in p.phases_into_bus[j]
                ks_on_phs = [k for k in j_to_k(j, p) if phs in p.phases_into_bus[k]]
                if !isempty(ks_on_phs)  # mid node
                    @constraint(m, [t in 1:p.Ntimesteps],
                        sum( Pᵢⱼ[string(i*"-"*j), phs, t] for i in i_to_j(j, p) ) +
                        Pⱼ[j,phs,t] - sum( Pᵢⱼ[string(j*"-"*k), phs, t] for k in ks_on_phs ) == 0
                    )
                    @constraint(m, [t in 1:p.Ntimesteps],
                        sum( Qᵢⱼ[string(i*"-"*j), phs, t] for i in i_to_j(j, p) ) +
                        Qⱼ[j,phs,t] - sum( Qᵢⱼ[string(j*"-"*k), phs, t] for k in ks_on_phs ) == 0
                    )
                else  # leaf node
                    @constraint(m, [t in 1:p.Ntimesteps],
                        sum( Pᵢⱼ[string(i*"-"*j), phs, t] for i in i_to_j(j, p) ) + Pⱼ[j,phs,t] == 0
                    )
                    @constraint(m, [t in 1:p.Ntimesteps],
                        sum( Qᵢⱼ[string(i*"-"*j), phs, t] for i in i_to_j(j, p) ) + Qⱼ[j,phs,t] == 0
                    )
                end

            end
        end
    end # intermediate and leaf nodes
end


function constrain_substation_voltage(m, p::Inputs)
    # @info "constrain_substation_voltage"
    @constraint(m, con_substationV[t in 1:p.Ntimesteps],
       m[:vsqrd][p.substation_bus, t] == p.v0^2
    )
    nothing
end


function constrain_substation_voltage(m, p::Inputs{ThreePhase})
    # @info "constrain_substation_voltage"
    @constraint(m, con_substationV[phs in 1:3, t in 1:p.Ntimesteps],
       m[:vsqrd][p.substation_bus, phs, t] == p.v0^2
    )
end


function constrain_KVL(m, p::Inputs)
    w = m[:vsqrd]
    P = m[:Pᵢⱼ]
    Q = m[:Qᵢⱼ]
    for j in p.busses
        for i in i_to_j(j, p)  # for radial network there is only one i in i_to_j
            i_j = string(i*"-"*j)
            rᵢⱼ = rij(i,j,p)
            xᵢⱼ = xij(i,j,p)
            vcon = @constraint(m, [t in 1:p.Ntimesteps],
                w[j,t] == w[i,t]
                    - 2*(rᵢⱼ * P[i_j,t] +  xᵢⱼ * Q[i_j,t])
            )
        end
    end
    nothing
end


function constrain_KVL(m, p::Inputs{ThreePhase})
    # @info "constrain_KVL"
    w = m[:vsqrd]
    P = m[:Pᵢⱼ]
    Q = m[:Qᵢⱼ]
    for j in p.busses
        for i in i_to_j(j, p)  # for radial network there is only one i in i_to_j
            i_j = string(i*"-"*j)
            MP = MPij(i,j,p)
            MQ = MQij(i,j,p)
            for phs in p.phases_into_bus[j]
                @constraint(m, [t in 1:p.Ntimesteps],
                    w[j,phs,t] == w[i,phs,t]
                        + sum(MP[phs,k] * P[i_j,k,t] for k=p.phases_into_bus[j])
                        + sum(MQ[phs,k] * Q[i_j,k,t] for k=p.phases_into_bus[j])
                )
            end
        end
    end
    # TODO need to set node w's to zero when a phase is not connected ?
end


"""
    constrain_loads(m, p::Inputs)

- set loads to negative of Inputs.Pload, which are normalized by Sbase when creating Inputs
- keys of Pload must match Inputs.busses. Any missing keys have load set to zero.
- Inputs.substation_bus is unconstrained, slack bus
"""
function constrain_loads(m, p::Inputs)
    Pⱼ = m[:Pⱼ]
    Qⱼ = m[:Qⱼ]
    
    for j in p.busses
        if j in keys(p.Pload)
            @constraint(m, [t in 1:p.Ntimesteps],
                Pⱼ[j,t] == -p.Pload[j][t] / p.Sbase
            )
        elseif j != p.substation_bus
            @constraint(m, [t in 1:p.Ntimesteps],
                Pⱼ[j,t] == 0
            )
        end
        if j in keys(p.Qload)
            @constraint(m, [t in 1:p.Ntimesteps],
                Qⱼ[j,t] == -p.Qload[j][t] / p.Sbase
            )
        elseif j != p.substation_bus
            @constraint(m, [t in 1:p.Ntimesteps],
                Qⱼ[j,t] == 0
            )
        end
    end
    nothing
end


"""
    constrain_loads(m, p::Inputs{ThreePhase})

- set loads to negative of Inputs.Pload, which are normalized by Sbase when creating Inputs
- keys of Pload must match Inputs.busses. Any missing keys have load set to zero.
- Inputs.substation_bus is unconstrained, slack bus
"""
function constrain_loads(m, p::Inputs{ThreePhase})
    Pⱼ = m[:Pⱼ]
    Qⱼ = m[:Qⱼ]
    a0 = 0.9  # TODO keep a0 and a1 values? (will have to update more P/Q constraints in 3 phase test)
    a1 = 0.1
    m[:cons] = Dict()
    m[:cons][:injection_equalities] = Dict()
    
    for j in p.busses
        m[:cons][:injection_equalities][j] = Dict()
        m[:cons][:injection_equalities][j][:P] = Dict()
        if j in keys(p.Pload)
            for phs in keys(p.Pload[j])
                if !(phs in p.phases_into_bus[j])
                    @warn "Load provided for bus $j, phase $phs but there are no lines into that point."
                else
                    m[:cons][:injection_equalities][j][:P][phs] = @constraint(m, [t in 1:p.Ntimesteps],
                        Pⱼ[j,phs,t] == -p.Pload[j][phs][t] / p.Sbase * (a0 + a1 * m[:vsqrd][j,phs,t])
                    )
                end
            end
            
            for phs in setdiff(p.phases_into_bus[j], keys(p.Pload[j]))
                m[:cons][:injection_equalities][j][:P][phs] = m[:cons][:injection_equalities][j][:P][phs] = @constraint(m, [t in 1:p.Ntimesteps],
                    Pⱼ[j,phs,t] == 0
                )
            end
        elseif j != p.substation_bus
            for phs in p.phases_into_bus[j]
                m[:cons][:injection_equalities][j][:P][phs] = @constraint(m, [t in 1:p.Ntimesteps],
                    Pⱼ[j,phs,t] == 0
                )
            end
        end  # real loads

        m[:cons][:injection_equalities][j][:Q] = Dict()
        if j in keys(p.Qload)
            for phs in keys(p.Qload[j])
                if !(phs in p.phases_into_bus[j])
                    @warn "Load provided for bus $j, phase $phs but there are no lines into that point."
                else
                    m[:cons][:injection_equalities][j][:Q][phs] = @constraint(m, [t in 1:p.Ntimesteps],
                        Qⱼ[j,phs,t] == -p.Qload[j][phs][t] / p.Sbase * (a0 + a1 * m[:vsqrd][j,phs,t])
                    )
                end
            end
            for phs in setdiff(p.phases_into_bus[j], keys(p.Qload[j]))
                m[:cons][:injection_equalities][j][:Q][phs] = @constraint(m, [t in 1:p.Ntimesteps],
                    Qⱼ[j,phs,t] == 0
                )
            end

        elseif j != p.substation_bus
            for phs in p.phases_into_bus[j]
                m[:cons][:injection_equalities][j][:Q][phs] = @constraint(m, [t in 1:p.Ntimesteps],
                    Qⱼ[j,phs,t] == 0
                )
            end
        end  # reactive loads
    end
end


function constrain_bounds(m::JuMP.AbstractModel, p::Inputs)
    @info("LinDistFlow.constrain_bounds is deprecated. Include bounds in LinDistFlow.Inputs.")
    nothing
end