# TODO make variable indices align with BranchFlowModel (and perhaps define some in CommonOPF?)
"""
    LinDistFlow.add_variables(m::JuMP.AbstractModel, p::Inputs{MultiPhase})

define: Pj, Qj, Pij, Qij, vsqrd

TODO:
- currently indexed bus/edge, phase, time
- but want to align with BranchFlowMmodel MultiPhase: time, bus/edge, phase (broadest to narrowest)
- and standardize variable containers in CommonOPF
- but using nested Dicts makes variable construction slow b/c have to create one at time,
    (which was necessary in MultiPhase BranchFlowModel) - but maybe not a big deal
"""
function add_variables(m::JuMP.AbstractModel, p::Inputs{MultiPhase})
    d = p.phases_into_bus
    d[p.substation_bus] = [1,2,3]
    T = 1:p.Ntimesteps
    # bus injections
    @variables m begin
        p.P_lo_bound <= Pj[b in p.busses, d[b], T] <= p.P_up_bound
        p.Q_lo_bound <= Qj[b in p.busses, d[b], T] <= p.Q_up_bound
    end
    
    # voltage squared
    @variable(m, p.v_lolim^2 <= vsqrd[b in p.busses, d[b], T] <= p.v_uplim^2 ) 
    
    edge2phases = Dict(k=>v for (k,v) in zip(p.edge_keys, p.phases))  # this will fail for mesh network

    # line flows, power sent from i to j
    @variable(m, p.P_lo_bound <= Pij[e in p.edge_keys, edge2phases[e], T] <= p.P_up_bound )
    
    @variable(m, p.Q_lo_bound <= Qij[e in p.edge_keys, edge2phases[e], T] <= p.Q_up_bound )
    # TODO line flow limit inputs
    nothing
end


function constrain_power_balance(m, p::Inputs{MultiPhase})
    Pj = m[:Pj]
    Qj = m[:Qj]
    Pij = m[:Pij]
    Qij = m[:Qij]

    source_nodes = setdiff([string(e[1]) for e in p.edges],[e[2] for e in p.edges])
    for j in source_nodes
        for phs in [1,2,3]
            ks_on_phs = [k for k in j_to_k(j, p) if phs in p.phases_into_bus[k]]
            if !isempty(ks_on_phs)
                @constraint(m, [t in 1:p.Ntimesteps],
                    Pj[j,phs,t] - sum( Pij[string(j*"-"*k), phs, t] for k in ks_on_phs ) == 0
                )
                @constraint(m, [t in 1:p.Ntimesteps],
                    Qj[j,phs,t] - sum( Qij[string(j*"-"*k), phs, t] for k in ks_on_phs ) == 0
                )
            else
                @warn "Source node $j is not connected to any busses on phase $phs. Setting Pj and Qj to zero."
                @constraint(m, [t in 1:p.Ntimesteps],
                    Pj[j,phs,t] == 0
                )
                @constraint(m, [t in 1:p.Ntimesteps],
                    Qj[j,phs,t] == 0
                )
            end
        end
    end # source nodes

    for j in setdiff(p.busses, source_nodes)
        if isempty(i_to_j(j, p)) && isempty(j_to_k(j, p))  # unconnected nodes
            @warn "Bus $j has no edges in or out; setting Pj and Qj to zero."
            @constraint(m, [phs in 1:3, t in 1:p.Ntimesteps],
                Pj[j,phs,t] == 0
            )
            @constraint(m, [phs in 1:3, t in 1:p.Ntimesteps],
                Qj[j,phs,t] == 0
            )
        else  
            for phs in p.phases_into_bus[j]
                ks_on_phs = [k for k in j_to_k(j, p) if phs in p.phases_into_bus[k]]
                if !isempty(ks_on_phs)  # mid node
                    @constraint(m, [t in 1:p.Ntimesteps],
                        sum( Pij[string(i*"-"*j), phs, t] for i in i_to_j(j, p) ) +
                        Pj[j,phs,t] - sum( Pij[string(j*"-"*k), phs, t] for k in ks_on_phs ) == 0
                    )
                    @constraint(m, [t in 1:p.Ntimesteps],
                        sum( Qij[string(i*"-"*j), phs, t] for i in i_to_j(j, p) ) +
                        Qj[j,phs,t] - sum( Qij[string(j*"-"*k), phs, t] for k in ks_on_phs ) == 0
                    )
                else  # leaf node
                    @constraint(m, [t in 1:p.Ntimesteps],
                        sum( Pij[string(i*"-"*j), phs, t] for i in i_to_j(j, p) ) + Pj[j,phs,t] == 0
                    )
                    @constraint(m, [t in 1:p.Ntimesteps],
                        sum( Qij[string(i*"-"*j), phs, t] for i in i_to_j(j, p) ) + Qj[j,phs,t] == 0
                    )
                end

            end
        end
    end # intermediate and leaf nodes
end


function constrain_substation_voltage(m, p::Inputs{MultiPhase})
    # @info "constrain_substation_voltage"
    @constraint(m, con_substationV[phs in 1:3, t in 1:p.Ntimesteps],
       m[:vsqrd][p.substation_bus, phs, t] == p.v0^2
    )
end


function constrain_KVL(m, p::Inputs{MultiPhase})
    # TODO store voltage constraints in model dict
    w = m[:vsqrd]
    P = m[:Pij]
    Q = m[:Qij]
    for j in p.busses
        for i in i_to_j(j, p)  # for radial network there is only one i in i_to_j
            if !( (i,j) in keys(p.regulators) )
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
            else
                if has_vreg(p, j)
                    for phs in p.phases_into_bus[j]
                        @constraint(m, [t in 1:p.Ntimesteps],
                            w[j,phs,t] == p.regulators[(i,j)][:vreg][phs]^2
                        )
                    end
                else  # default turn_ratio is 1.0
                    for phs in p.phases_into_bus[j]
                        @constraint(m, [t in 1:p.Ntimesteps],
                            w[j,phs,t] == w[i,phs,t] * p.regulators[(i,j)][:turn_ratio][phs]^2 
                        )
                    end
                end
            end
        end
    end
    # TODO need to set node w's to zero when a phase is not connected ?
end


"""
    constrain_loads(m, p::Inputs{MultiPhase})

- set loads to negative of Inputs.Pload and Qload, which are normalized by Sbase when creating Inputs
- keys of Pload and Qload must match Inputs.busses. Any missing keys have load set to zero.
- Inputs.substation_bus is unconstrained, slack bus
"""
function constrain_loads(m, p::Inputs{MultiPhase})
    Pj = m[:Pj]
    Qj = m[:Qj]
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
                        Pj[j,phs,t] == -p.Pload[j][phs][t] / p.Sbase
                    )
                end
            end
            
            for phs in setdiff(p.phases_into_bus[j], keys(p.Pload[j]))
                m[:cons][:injection_equalities][j][:P][phs] = m[:cons][:injection_equalities][j][:P][phs] = @constraint(m, [t in 1:p.Ntimesteps],
                    Pj[j,phs,t] == 0
                )
            end
        elseif j != p.substation_bus
            for phs in p.phases_into_bus[j]
                m[:cons][:injection_equalities][j][:P][phs] = @constraint(m, [t in 1:p.Ntimesteps],
                    Pj[j,phs,t] == 0
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
                        Qj[j,phs,t] == -p.Qload[j][phs][t] / p.Sbase
                    )
                end
            end
            for phs in setdiff(p.phases_into_bus[j], keys(p.Qload[j]))
                m[:cons][:injection_equalities][j][:Q][phs] = @constraint(m, [t in 1:p.Ntimesteps],
                    Qj[j,phs,t] == 0
                )
            end

        elseif j != p.substation_bus
            for phs in p.phases_into_bus[j]
                m[:cons][:injection_equalities][j][:Q][phs] = @constraint(m, [t in 1:p.Ntimesteps],
                    Qj[j,phs,t] == 0
                )
            end
        end  # reactive loads
    end
end


function define_line_amp_estimates(m::JuMP.AbstractModel, p::Inputs{MultiPhase})
    Pij = m[:Pij]
    Qij = m[:Qij]
    m[:amps_pu] = Dict(ek => Dict() for ek in p.edge_keys)
    for j in p.busses
        for i in i_to_j(j, p)  # for radial network there is only one i in i_to_j
            i_j = string(i*"-"*j)
            r = rij(i,j,p)
            x = xij(i,j,p)
            z = sqrt(r^2 + x^2)  # in per-unit
            for phs in p.phases_into_bus[j]
                m[:amps_pu][i*"-"*j][phs] = @expression(m, [t in 1:p.Ntimesteps],
                    sum(r[phs,k] * Pij[i_j,k,t] + x[phs,k] * Qij[i_j,k,t] for k=p.phases_into_bus[j]) / z[phs,phs]
                )
            end
        end
    end
end


"""
    constrain_line_amps(m::JuMP.AbstractModel, p::Inputs{MultiPhase})

Estimating the line amps as ``|(V_i - V_j) / Z|`` where we use the approximation:

``
|V_i - V_j| \\approx r_{ij} P_{ij} + x_{ij} Q_{ij}
``

The expressions are stored in the `model.obj_dict` with key `:amps_pu`.
"""
function constrain_line_amps(m::JuMP.AbstractModel, p::Inputs{MultiPhase})
    if !(:amps_pu in keys(m.obj_dict))
        define_line_amp_estimates(m, p)
    end
    for j in p.busses
        for i in i_to_j(j, p)  # for radial network there is only one i in i_to_j
            ij_linecode = get_ijlinecode(i,j,p)
            amps_pu_limit = sqrt(p.Isquared_up_bounds[ij_linecode]) / p.Ibase
            for phs in p.phases_into_bus[j]
                @constraint(m, [t in 1:p.Ntimesteps],
                    -amps_pu_limit <= m[:amps_pu][i*"-"*j][phs][t] <= amps_pu_limit
                )
            end
        end
    end
end
