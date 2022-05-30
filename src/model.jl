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

    p.Nlte_cons += 6 * p.Nnodes * p.Ntimesteps
    
    ij_edges = [string(i*"-"*j) for j in p.busses for i in i_to_j(j, p)]
    Nedges = length(ij_edges)

    # line flows, power sent from i to j
    @variable(m, p.P_lo_bound <= Pᵢⱼ[ij_edges, T] <= p.P_up_bound )
    
    @variable(m, p.Q_lo_bound <= Qᵢⱼ[ij_edges, T] <= p.Q_up_bound )

    p.Nlte_cons += 4 * Nedges * p.Ntimesteps
end


function add_variables(m, p::Inputs{ThreePhase})
    T = 1:p.Ntimesteps
    # bus injections
    @variables m begin
        p.P_lo_bound <= Pⱼ[p.busses, 1:3, T] <= p.P_up_bound
        p.Q_lo_bound <= Qⱼ[p.busses, 1:3, T] <= p.Q_up_bound
    end

    # voltage squared
    @variable(m, p.v_lolim^2 <= vsqrd[p.busses, 1:3, T] <= p.v_uplim^2 ) 

    p.Nlte_cons += 6 * p.Nnodes * p.Ntimesteps * 3
    
    ij_edges = [string(i*"-"*j) for j in p.busses for i in i_to_j(j, p)]
    Nedges = length(ij_edges)

    # line flows, power sent from i to j
    @variable(m, p.P_lo_bound <= Pᵢⱼ[ij_edges, 1:3, T] <= p.P_up_bound )
    
    @variable(m, p.Q_lo_bound <= Qᵢⱼ[ij_edges, 1:3, T] <= p.Q_up_bound )

    p.Nlte_cons += 4 * Nedges * p.Ntimesteps * 3
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
    p.Nequality_cons += 2 * p.Nnodes * p.Ntimesteps
end


function constrain_power_balance(m, p::Inputs{ThreePhase})
    Pⱼ = m[:Pⱼ]
    Qⱼ = m[:Qⱼ]
    Pᵢⱼ = m[:Pᵢⱼ]
    Qᵢⱼ = m[:Qᵢⱼ]
    # TODO change Pⱼ and Qⱼ to expressions, make P₀ and Q₀ dv's, which will reduce # of variables
    # by (Nnodes - 1)*8760 and number of constraints by 6*(Nnodes - 1)*8760
    for j in p.busses
        if isempty(i_to_j(j, p)) && !isempty(j_to_k(j, p)) # source nodes
            pcon = @constraint(m,  [phs in 1:3, t in 1:p.Ntimesteps],
                Pⱼ[j,phs,t] - sum( Pᵢⱼ[string(j*"-"*k), phs, t] for k in j_to_k(j, p) ) == 0
            )
            qcon = @constraint(m, [phs in 1:3, t in 1:p.Ntimesteps],
                Qⱼ[j,phs,t] - sum( Qᵢⱼ[string(j*"-"*k), phs, t] for k in j_to_k(j, p) ) == 0
            )
        elseif isempty(i_to_j(j, p)) && isempty(j_to_k(j, p))  # unconnected nodes
            @warn "Bus $j has no edges, setting Pⱼ and Qⱼ to zero."
            pcon = @constraint(m, [t in 1:p.Ntimesteps],
                Pⱼ[j,phs,t] == 0
            )
            qcon = @constraint(m, [phs in 1:3, t in 1:p.Ntimesteps],
                Qⱼ[j,phs,t] == 0
            )
        elseif !isempty(i_to_j(j, p)) && isempty(j_to_k(j, p))  # leaf nodes
            pcon = @constraint(m, [phs in 1:3, t in 1:p.Ntimesteps],
                sum( Pᵢⱼ[string(i*"-"*j), phs, t] for i in i_to_j(j, p) ) + Pⱼ[j, phs, t] == 0
            )
            qcon = @constraint(m, [phs in 1:3, t in 1:p.Ntimesteps],
                sum( Qᵢⱼ[string(i*"-"*j), phs, t] for i in i_to_j(j, p) ) + Qⱼ[j, phs, t] == 0
            )
        else
            pcon =  @constraint(m, [phs in 1:3, t in 1:p.Ntimesteps],
                sum( Pᵢⱼ[string(i*"-"*j), phs, t] for i in i_to_j(j, p) ) +
                Pⱼ[j,phs,t] - sum( Pᵢⱼ[string(j*"-"*k), phs, t] for k in j_to_k(j, p) ) == 0
            )
            qcon = @constraint(m, [phs in 1:3, t in 1:p.Ntimesteps],
                sum( Qᵢⱼ[string(i*"-"*j), phs, t] for i in i_to_j(j, p) ) +
                Qⱼ[j,phs,t] - sum( Qᵢⱼ[string(j*"-"*k), phs, t] for k in j_to_k(j, p) ) == 0
            )
        end
    end
    p.Nequality_cons += 2 * p.Nnodes * p.Ntimesteps *3
end


function constrain_substation_voltage(m, p::Inputs)
    # @info "constrain_substation_voltage"
    @constraint(m, con_substationV[t in 1:p.Ntimesteps],
       m[:vsqrd][p.substation_bus, t] == p.v0^2
    )
    p.Nequality_cons += p.Ntimesteps
end


function constrain_substation_voltage(m, p::Inputs{ThreePhase})
    # @info "constrain_substation_voltage"
    @constraint(m, con_substationV[phs in 1:3, t in 1:p.Ntimesteps],
       m[:vsqrd][p.substation_bus, phs, t] == p.v0^2
    )
    p.Nequality_cons += p.Ntimesteps * 3
end


function constrain_KVL(m, p::Inputs)
    # @info "constrain_KVL"
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
    p.Nequality_cons += length(p.edges) * p.Ntimesteps
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
            @constraint(m, [phs in 1:3, t in 1:p.Ntimesteps],
                w[j,phs,t] == w[i,phs,t]
                    - sum(MP[phs,k] * P[i_j,k,t] for k=1:3)
                    - sum(MQ[phs,k] * Q[i_j,k,t] for k=1:3)
            )
        end
    end
    # TODO need to set node w's to zero when a phase is not connected ?
    p.Nequality_cons += length(p.edges) * p.Ntimesteps * 3
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
    
    p.Nequality_cons += 2 * (p.Nnodes - 1) * p.Ntimesteps
end


"""
    constrain_loads(m, p::Inputs)

- set loads to negative of Inputs.Pload, which are normalized by Sbase when creating Inputs
- keys of Pload must match Inputs.busses. Any missing keys have load set to zero.
- Inputs.substation_bus is unconstrained, slack bus
"""
function constrain_loads(m, p::Inputs{ThreePhase})
    Pⱼ = m[:Pⱼ]
    Qⱼ = m[:Qⱼ]
    
    for j in p.busses
        if j in keys(p.Pload)
            for phs in keys(p.Pload[j])
                @constraint(m, [t in 1:p.Ntimesteps],
                    Pⱼ[j,phs,t] == -p.Pload[j][phs][t] / p.Sbase
                )
            end
            # for phs in setdiff([1,2,3], keys(p.Pload[j]))
            #     @constraint(m, [t in 1:p.Ntimesteps],
            #         Pⱼ[j,phs,t] == 0
            #     )
            # end
        elseif j != p.substation_bus
            @constraint(m, [phs in 1:3, t in 1:p.Ntimesteps],
                Pⱼ[j,phs,t] == 0
            )
        end
        if j in keys(p.Qload)
            for phs in keys(p.Qload[j])
                @constraint(m, [t in 1:p.Ntimesteps],
                    Qⱼ[j,phs,t] == -p.Qload[j][phs][t] / p.Sbase
                )
            end
            # for phs in setdiff([1,2,3], keys(p.Qload[j]))
            #     @constraint(m, [t in 1:p.Ntimesteps],
            #         Qⱼ[j,phs,t] == 0
            #     )
            # end
        elseif j != p.substation_bus
            @constraint(m, [phs in 1:3, t in 1:p.Ntimesteps],
                Qⱼ[j,phs,t] == 0
            )
        end
    end
    
    p.Nequality_cons += 2 * (p.Nnodes - 1) * p.Ntimesteps * 3
end


function constrain_bounds(m::JuMP.AbstractModel, p::Inputs)
    @info("LinDistFlow.constrain_bounds is deprecated. Include bounds in LinDistFlow.Inputs.")
    nothing
end