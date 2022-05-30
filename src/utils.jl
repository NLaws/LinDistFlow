"""
    function i_to_j(j::String, p::Inputs)
find all busses upstream of bus j

!!! note
    In a radial network this function should return an Array with length of 1.
"""
function i_to_j(j::String, p::Inputs)
    convert(Array{String, 1}, map(x->x[1], filter(t->t[2]==j, p.edges)))
end


"""
    function j_to_k(j::String, p::Inputs)
find all busses downstream of bus j
"""
function j_to_k(j::String, p::Inputs)
    convert(Array{String, 1}, map(x->x[2], filter(t->t[1]==j, p.edges)))
end


function rij(i::String, j::String, p::Inputs{SinglePhase})
    linecode = get_ijlinecode(i, j, p)
    linelength = get_ijlinelength(i, j, p)
    rmatrix = p.Zdict[linecode]["rmatrix"] * linelength / p.Zbase
    return rmatrix[1]
end


function rij(i::String, j::String, p::Inputs{ThreePhase})
    linecode = get_ijlinecode(i, j, p)
    linelength = get_ijlinelength(i, j, p)
    rmatrix = p.Zdict[linecode]["rmatrix"] * linelength / p.Zbase
    return rmatrix
end


function xij(i::String, j::String, p::Inputs{SinglePhase})
    linecode = get_ijlinecode(i, j, p)
    linelength = get_ijlinelength(i, j, p)
    xmatrix = p.Zdict[linecode]["xmatrix"] * linelength / p.Zbase
    return xmatrix[1]
end


function xij(i::String, j::String, p::Inputs{ThreePhase})
    linecode = get_ijlinecode(i, j, p)
    linelength = get_ijlinelength(i, j, p)
    xmatrix = p.Zdict[linecode]["xmatrix"] * linelength / p.Zbase
    return xmatrix
end


function get_ijlinelength(i::String, j::String, p::Inputs)
    ij_idxs = get_ij_idxs(i, j, p)
    return p.linelengths[ij_idxs[1]]
end


function get_ijlinecode(i::String, j::String, p::Inputs)
    ij_idxs = get_ij_idxs(i, j, p)
    return p.linecodes[ij_idxs[1]]
end


function get_ijedge(i::String, j::String, p::Inputs)
    ij_idxs = get_ij_idxs(i, j, p)
    return p.edges[ij_idxs[1]]
end


function get_ij_idxs(i::String, j::String, p::Inputs)
    ij_idxs = findall(t->(t[1]==i && t[2]==j), p.edges)
    if length(ij_idxs) > 1
        error("found more than one edge for i=$i and j=$j")
    elseif length(ij_idxs) == 0
        error("found no matching edges for i=$i and j=$j")
    else
        return ij_idxs
    end
end


function get_edge_values(var_prefix::String, m::JuMP.AbstractModel, p::Inputs)
    vals = Float64[]
    for edge in p.edges
        var = string(var_prefix, "[", edge[1], "-", edge[2], "]")
        try
            val = value(variable_by_name(m, var))
            if startswith(var_prefix, "l")
                val = sqrt(val)
            end
            push!(vals, round(val; digits=8))
        catch e
            println(var, "failed", e)
        end
    end
    return vals
end


function get_bus_values(var_prefix::String, m::JuMP.AbstractModel, p::Inputs)
    vals = Float64[]
    for b in p.busses
        var = string(var_prefix,  "[", b, "]")
        try
            val = value(variable_by_name(m, var))
            if startswith(var_prefix, "v")
                val = sqrt(val)
            end
            push!(vals, round(val; digits=7))
        catch e
            println(var, " failed: ", e)
        end
    end
    return vals
end


function get_constraints_by_variable_name(m, v::String)
    ac = ConstraintRef[]
    for tup in list_of_constraint_types(m)
        append!(ac, all_constraints(m, tup[1], tup[2]))
    end
    filter( cr -> occursin(v, string(cr)), ac )
end


"""

Real power coefficients for 3 phase voltage drop from node i to j
"""
function MPij(i::String, j::String, p::Inputs{ThreePhase})
    M = zeros((3,3))
    r = rij(i,j,p)
    x = xij(i,j,p)
    return M
end


"""
    function recover_voltage_current(m, p::Inputs)

Algorithm 2 from Gan & Low 2014
"""
function recover_voltage_current(m, p::Inputs, nodetobusphase)
    # TODO finish converting this to single phase, validate
    Iij = Dict{String, Array{Complex, 1}}()
    Vj = Dict{String, Array{Complex, 1}}()
    Vj[p.substation_bus] = [p.v0*exp(0im)]

    for i in p.busses
        for j in j_to_k(i,p)
            tr_vᵢ = zeros(3)
            Pᵢⱼ = zeros(3,3)
            Qᵢⱼ = zeros(3,3)

            tr_vᵢ = value(variable_by_name(m, string("vⱼ", "[", i, "]")))
            Pᵢⱼ = value(variable_by_name(m, string("Pᵢⱼ", "[", i, "-", j, "]")))
            Qᵢⱼ = value(variable_by_name(m, string("Qᵢⱼ", "[", i, "-", j, "]")))
            
            Sᵢⱼ = @. complex(Pᵢⱼ, Qᵢⱼ)
            r = rij(i,j,p)
            x = xij(i,j,p)
            zᵢⱼ = @. complex(r, x)

            Iij[i*"-"*j] = 1/sum(tr_vᵢ) * Sᵢⱼ' * Vj[i]
            Vj[j] = Vj[i] - zᵢⱼ * Iij[i*"-"*j]
        end
    end
    # convert dicts to vectors in same order as nodetobusphase
    v = Complex[]
    c = Complex[]

    for busphase in [split(bp, ".") for bp in nodetobusphase]
        b = busphase[1]
        ph = parse(Int, busphase[2])
        push!(v, Vj[b][ph])
    end
    return v, Iij
end

# TODO add check of rank 1 for exact solution (10g) in Gan and Low