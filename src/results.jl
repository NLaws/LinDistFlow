abstract type AbstractResults end

"""
TODO
"""
struct Results <: AbstractResults
    voltage_magnitudes
    real_power_injections
    reactive_power_injections
    real_sending_end_powers
    reactive_sending_end_powers
    shadow_prices
end


"""
    Results(m::AbstractModel, p::Inputs{SinglePhase})

return a `Results` struct with fieldnames:

    voltage_magnitudes
    real_power_injections
    reactive_power_injections
    current_magnitudes
    real_sending_end_powers
    reactive_sending_end_powers
    prices

"""
function Results(m::AbstractModel, p::Inputs{SinglePhase}; digits=8)

    vs  = CommonOPF.get_variable_values(:vsqrd, m, p, digits=digits)
    Pj  = CommonOPF.get_variable_values(:Pj,    m, p, digits=digits)
    Qj  = CommonOPF.get_variable_values(:Qj,    m, p, digits=digits)
    Pij = CommonOPF.get_variable_values(:Pij,   m, p, digits=digits)
    Qij = CommonOPF.get_variable_values(:Qij,   m, p, digits=digits)
    prices = Dict(
        j => JuMP.dual.(m[:loadbalcons][j]["p"])
        for j in p.busses
    )

    Results(vs, Pj, Qj, Pij, Qij, prices)
end


"""
    Results(m::AbstractModel, p::Inputs{MultiPhase})

return a `Results` struct with fieldnames:

    voltage_magnitudes
    real_power_injections
    reactive_power_injections
    current_magnitudes
    real_sending_end_powers
    reactive_sending_end_powers

"""
function Results(m::AbstractModel, p::Inputs{MultiPhase}; digits=8)

    vs  = CommonOPF.get_variable_values(:vsqrd, m, p, digits=digits)
    Pj  = CommonOPF.get_variable_values(:Pj,    m, p, digits=digits)
    Qj  = CommonOPF.get_variable_values(:Qj,    m, p, digits=digits)
    Pij = CommonOPF.get_variable_values(:Pij,   m, p, digits=digits)
    Qij = CommonOPF.get_variable_values(:Qij,   m, p, digits=digits)
    prices = Dict("TODO" => "for MultiPhase")

    Results(vs, Pj, Qj, Pij, Qij, prices)
end


"""
    get_variable_values(var::Symbol, m::JuMP.AbstractModel, p::Inputs{MultiPhase}; digits=6)

!!! note
    This method only works with LinDistFlow MultiPhase models (not BranchFlowModel MultiPhase)

!!! note
    Rounding can be necessary for values that require `sqrt` and have optimal values of zero like 
    `-3.753107219618953e-31`
"""
function CommonOPF.get_variable_values(var::Symbol, m::JuMP.AbstractModel, p::Inputs{MultiPhase}; digits=8)
    d = Dict()
    if var in [:Pj, :Qj, :vsqrd]  # TODO make these a const in CommonOPF
        vals = value.(m[var])
        for b in p.busses
            d[b] = Dict()
            for phs in p.phases_into_bus[b]
                d[b][phs] = round.(values(vals[b,phs,:].data), digits=digits)
                if var == :vsqrd
                    d[b][phs] = sqrt.(d[b][phs])
                else
                    d[b][phs] *= p.Sbase  # scale powers back to absolute units TODO in BFM
                end
            end
        end
    elseif var in [:Pij, :Qij, :lij]  # TODO make these a const in CommonOPF TODO in BFM
        vals = value.(m[var])
        for (edge_key, phases) in zip(p.edge_keys, p.phases)
            d[edge_key] = Dict()
            for phs in phases
                d[edge_key][phs] = round.(values(vals[edge_key,phs,:].data), digits=digits)
                if var == :lij
                    d[edge_key][phs] = sqrt.(d[edge_key][phs])
                else
                    d[edge_key][phs] *= p.Sbase  # scale powers back to absolute units
                end
            end
        end
    else
        @warn "$var is not a valid variable symbol"
    end
    return d
end