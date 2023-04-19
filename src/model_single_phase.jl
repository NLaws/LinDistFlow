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
    # TODO store constraints in m.obj_dict in standardized fashion
end


function add_variables(m, p::Inputs)
    T = 1:p.Ntimesteps
    # bus injections
    @variables m begin
        p.P_lo_bound <= Pj[p.busses, T] <= p.P_up_bound
        p.Q_lo_bound <= Qj[p.busses, T] <= p.Q_up_bound
    end

    # voltage squared
    @variable(m, p.v_lolim^2 <= vsqrd[p.busses, T] <= p.v_uplim^2 ) 
    
    # line flows, power sent from i to j
    @variable(m, p.P_lo_bound <= Pij[p.edge_keys, T] <= p.P_up_bound )
    
    @variable(m, p.Q_lo_bound <= Qij[p.edge_keys, T] <= p.Q_up_bound )

    nothing
end


function constrain_power_balance(m, p::Inputs)
    Pj = m[:Pj]
    Qj = m[:Qj]
    Pij = m[:Pij]
    Qij = m[:Qij]
    # TODO change Pj and Qj to expressions, make P₀ and Q₀ dv's, which will reduce # of variables
    # by (Nnodes - 1)*8760 and number of constraints by 6*(Nnodes - 1)*8760
    for j in p.busses
        if isempty(i_to_j(j, p)) && !isempty(j_to_k(j, p)) # source nodes
            pcon = @constraint(m,  [t in 1:p.Ntimesteps],
                Pj[j,t] - sum( Pij[string(j*"-"*k), t] for k in j_to_k(j, p) ) == 0
            )
            qcon = @constraint(m, [t in 1:p.Ntimesteps],
                Qj[j,t] - sum( Qij[string(j*"-"*k), t] for k in j_to_k(j, p) ) == 0
            )
        elseif isempty(i_to_j(j, p)) && isempty(j_to_k(j, p))  # unconnected nodes
            @warn "Bus $j has no edges, setting Pj and Qj to zero."
            pcon = @constraint(m, [t in 1:p.Ntimesteps],
                Pj[j,t] == 0
            )
            qcon = @constraint(m, [t in 1:p.Ntimesteps],
                Qj[j,t] == 0
            )
        elseif !isempty(i_to_j(j, p)) && isempty(j_to_k(j, p))  # leaf nodes
            pcon = @constraint(m, [t in 1:p.Ntimesteps],
                sum( Pij[string(i*"-"*j), t] for i in i_to_j(j, p) ) + Pj[j, t] == 0
            )
            qcon = @constraint(m, [t in 1:p.Ntimesteps],
                sum( Qij[string(i*"-"*j), t] for i in i_to_j(j, p) ) + Qj[j, t] == 0
            )
        else
            pcon =  @constraint(m, [t in 1:p.Ntimesteps],
                sum( Pij[string(i*"-"*j), t] for i in i_to_j(j, p) ) +
                Pj[j,t] - sum( Pij[string(j*"-"*k), t] for k in j_to_k(j, p) ) == 0
            )
            qcon = @constraint(m, [t in 1:p.Ntimesteps],
                sum( Qij[string(i*"-"*j), t] for i in i_to_j(j, p) ) +
                Qj[j,t] - sum( Qij[string(j*"-"*k), t] for k in j_to_k(j, p) ) == 0
            )
        end
    end
    nothing
end


function constrain_substation_voltage(m, p::Inputs)
    # @info "constrain_substation_voltage"
    @constraint(m, con_substationV[t in 1:p.Ntimesteps],
       m[:vsqrd][p.substation_bus, t] == p.v0^2
    )
    nothing
end


function constrain_KVL(m, p::Inputs)
    w = m[:vsqrd]
    P = m[:Pij]
    Q = m[:Qij]
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


"""
    constrain_loads(m, p::Inputs)

- set loads to negative of Inputs.Pload, which are normalized by Sbase when creating Inputs
- keys of Pload must match Inputs.busses. Any missing keys have load set to zero.
- Inputs.substation_bus is unconstrained, slack bus
"""
function constrain_loads(m, p::Inputs)
    Pj = m[:Pj]
    Qj = m[:Qj]
    
    for j in p.busses
        if j in keys(p.Pload)
            @constraint(m, [t in 1:p.Ntimesteps],
                Pj[j,t] == -p.Pload[j][t] / p.Sbase
            )
        elseif j != p.substation_bus
            @constraint(m, [t in 1:p.Ntimesteps],
                Pj[j,t] == 0
            )
        end
        if j in keys(p.Qload)
            @constraint(m, [t in 1:p.Ntimesteps],
                Qj[j,t] == -p.Qload[j][t] / p.Sbase
            )
        elseif j != p.substation_bus
            @constraint(m, [t in 1:p.Ntimesteps],
                Qj[j,t] == 0
            )
        end
    end
    nothing
end
