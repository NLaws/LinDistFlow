function build_ldf!(m::JuMP.AbstractModel, p::Inputs)

    add_variables(m, p)
    constrain_power_balance(m, p)  # (10a)
    constrain_substation_voltage(m, p)  # (10c)
    constrain_KVL(m, p)  # (10e)
    constrain_loads(m, p)
    constrain_bounds(m, p)

end


function add_variables(m, p::Inputs)
    T = 1:p.Ntimesteps
    # bus injections
    @variables m begin
        Pⱼ[p.busses, T]
        Qⱼ[p.busses, T]
    end

    # voltage squared
    @variable(m, 0 <= vsqrd[p.busses, T])
    
    ij_edges = [string(i*"-"*j) for j in p.busses for i in i_to_j(j, p)]

    # line flows, power sent from i to j
    @variable(m, Pᵢⱼ[ij_edges, T])
    
    @variable(m, Qᵢⱼ[ij_edges, T])
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


function constrain_substation_voltage(m, p::Inputs)
    # @info "constrain_substation_voltage"
    @constraint(m, con_substationV[t in 1:p.Ntimesteps],
       m[:vsqrd][p.substation_bus, t] == p.v0^2
    )
    p.Nequality_cons += p.Ntimesteps
end


function constrain_KVL(m, p::Inputs)
    # @info "constrain_KVL"
    w = m[:vsqrd]
    P = m[:Pᵢⱼ]
    Q = m[:Qᵢⱼ]
    for j in p.busses
        for i in i_to_j(j, p)
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
                Pⱼ[j,t] == -p.Pload[j][t]
            )
        elseif j != p.substation_bus
            @constraint(m, [t in 1:p.Ntimesteps],
                Pⱼ[j,t] == 0
            )
        end
        if j in keys(p.Qload)
            @constraint(m, [t in 1:p.Ntimesteps],
                Qⱼ[j,t] == -p.Qload[j][t]
            )
        elseif j != p.substation_bus
            @constraint(m, [t in 1:p.Ntimesteps],
                Qⱼ[j,t] == 0
            )
        end
    end
    
    p.Nequality_cons += 2 * (p.Nnodes - 1) * p.Ntimesteps
end


function constrain_bounds(m::JuMP.AbstractModel, p::Inputs)
    # TODO add bounds to Inputs
    # ld_busses = setdiff(p.busses, [p.substation_bus])
    T = 1:p.Ntimesteps
    
    # node injections
    @constraint(m, [b in p.busses, t in T],
        m[:Pⱼ][b, t] <= 1e4
    )
    p.Nlte_cons += p.Nnodes * p.Ntimesteps
    
    @constraint(m, [b in p.busses, t in T],
        -1e4 - m[:Pⱼ][b, t] <= 0
    )
    p.Nlte_cons += p.Nnodes * p.Ntimesteps
    
    @constraint(m, [b in p.busses, t in T],
        m[:Qⱼ][b, t] <= 1e4
    )
    p.Nlte_cons += p.Nnodes * p.Ntimesteps
    
    @constraint(m, [b in p.busses, t in T],
        -1e4 - m[:Qⱼ][b, t] <= 0
    )
    p.Nlte_cons += p.Nnodes * p.Ntimesteps
    
    # voltage squared
    @constraint(m, [b in p.busses, t in T],
        m[:vsqrd][b, t] <= p.v_uplim^2
    )
    p.Nlte_cons += p.Nnodes * p.Ntimesteps
    
    @constraint(m, [b in p.busses, t in T],
        p.v_lolim^2 - m[:vsqrd][b, t] <= 0
    )
    p.Nlte_cons += p.Nnodes * p.Ntimesteps
    
    # line flows, power sent from i to j
    ij_edges = [string(i*"-"*j) for j in p.busses for i in i_to_j(j, p)]
    Nedges = length(ij_edges)
    
    @constraint(m, [e in ij_edges, t in T],
        m[:Pᵢⱼ][e,t] ≤ 1e4
    )
    p.Nlte_cons += Nedges * p.Ntimesteps
    
    @constraint(m, [e in ij_edges, t in T],
        -1e4 - m[:Pᵢⱼ][e,t] ≤ 0
    )
    p.Nlte_cons += Nedges * p.Ntimesteps
    
    @constraint(m, [e in ij_edges, t in T],
        m[:Qᵢⱼ][e,t] ≤ 1e4
    )
    p.Nlte_cons += Nedges * p.Ntimesteps
    
    @constraint(m, [e in ij_edges, t in T],
        -1e4 - m[:Qᵢⱼ][e,t] ≤ 0
    )
    p.Nlte_cons += Nedges * p.Ntimesteps
    nothing
end