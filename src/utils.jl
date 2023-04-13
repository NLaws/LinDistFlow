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


function rij(i::String, j::String, p::Inputs{MultiPhase})
    linecode = get_ijlinecode(i, j, p)
    linelength = get_ijlinelength(i, j, p)
    raw_rmatrix = p.Zdict[linecode]["rmatrix"] * linelength / p.Zbase
    phases = get_ijphases(i, j, p)
    if length(phases) == 3
        return raw_rmatrix
    end
    rmatrix = zeros(3,3)
    if length(phases) == 1
        rmatrix[phases[1], phases[1]] = raw_rmatrix[1]
        return rmatrix
    end
    # else have to handle two phases, which could be [1,2], [1,3], or [2,3]
    phs1, phs2 = phases
    rmatrix[phs1,phs1] = raw_rmatrix[1,1]
    rmatrix[phs2,phs2] = raw_rmatrix[2,2]
    rmatrix[phs1,phs2] = raw_rmatrix[1,2]
    rmatrix[phs2,phs1] = raw_rmatrix[2,1]
    return rmatrix
end


function xij(i::String, j::String, p::Inputs{SinglePhase})
    linecode = get_ijlinecode(i, j, p)
    linelength = get_ijlinelength(i, j, p)
    xmatrix = p.Zdict[linecode]["xmatrix"] * linelength / p.Zbase
    return xmatrix[1]
end


function xij(i::String, j::String, p::Inputs{MultiPhase})
    linecode = get_ijlinecode(i, j, p)
    linelength = get_ijlinelength(i, j, p)
    raw_xmatrix = p.Zdict[linecode]["xmatrix"] * linelength / p.Zbase
    phases = get_ijphases(i, j, p)
    if length(phases) == 3
        return raw_xmatrix
    end
    xmatrix = zeros(3,3)
    if length(phases) == 1
        xmatrix[phases[1], phases[1]] = raw_xmatrix[1]
        return xmatrix
    end
    # else have to handle two phases, which could be [1,2], [1,3], or [2,3]
    phs1, phs2 = phases
    xmatrix[phs1,phs1] = raw_xmatrix[1,1]
    xmatrix[phs2,phs2] = raw_xmatrix[2,2]
    xmatrix[phs1,phs2] = raw_xmatrix[1,2]
    xmatrix[phs2,phs1] = raw_xmatrix[2,1]
    return xmatrix
end


function get_ijlinelength(i::String, j::String, p::Inputs)
    ij_idx = get_ij_idx(i, j, p)
    return p.linelengths[ij_idx]
end


function get_ijlinecode(i::String, j::String, p::Inputs)
    ij_idx = get_ij_idx(i, j, p)
    return p.linecodes[ij_idx]
end


function get_ijedge(i::String, j::String, p::Inputs)
    ij_idx = get_ij_idx(i, j, p)
    return p.edges[ij_idx]
end


function get_ijphases(i::String, j::String, p::Inputs{MultiPhase})
    ij_idx = get_ij_idx(i, j, p)
    return p.phases[ij_idx]
end


function get_ij_idx(i::String, j::String, p::Inputs)
    ij_idxs = findall(t->(t[1]==i && t[2]==j), p.edges)
    if length(ij_idxs) > 1
        error("found more than one edge for i=$i and j=$j")
    elseif length(ij_idxs) == 0
        error("found no matching edges for i=$i and j=$j")
    else
        return ij_idxs[1]
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


function get_bus_values(var_prefix::String, m::JuMP.AbstractModel, p::Inputs{MultiPhase})
    var_refs = collect(v for v in all_variables(m) if startswith(string(v), var_prefix))
    d = Dict{String, Real}()
    for vr in var_refs
        d[string(vr)] = value(vr)
    end
    return d
end


function get_constraints_by_variable_name(m, v::String)
    ac = ConstraintRef[]
    for tup in list_of_constraint_types(m)
        append!(ac, all_constraints(m, tup[1], tup[2]))
    end
    filter( cr -> occursin(v, string(cr)), ac )
end


"""
    MPij(i::String, j::String, p::Inputs{MultiPhase})

Real power coefficients for 3 phase voltage drop from node i to j
"""
function MPij(i::String, j::String, p::Inputs{MultiPhase})
    M = zeros((3,3))
    r = rij(i,j,p)
    x = xij(i,j,p)
    M[1,:] = [-2r[1,1]             r[1,2]-sqrt(3)x[1,2] r[1,3]+sqrt(3)x[1,3]]
    M[2,:] = [r[2,1]+sqrt(3)x[2,1] -2r[2,2]             r[2,3]-sqrt(3)x[2,3]]
    M[3,:] = [r[3,1]-sqrt(3)x[3,1] r[3,2]+sqrt(3)x[3,2] -2r[3,3]            ]
    return M
end


"""
    MQij(i::String, j::String, p::Inputs{MultiPhase})

Reactive power coefficients for 3 phase voltage drop from node i to j
"""
function MQij(i::String, j::String, p::Inputs{MultiPhase})
    M = zeros((3,3))
    r = rij(i,j,p)
    x = xij(i,j,p)
    M[1,:] = [-2x[1,1]             x[1,2]+sqrt(3)r[1,2] x[1,3]-sqrt(3)r[1,3]]
    M[2,:] = [x[2,1]-sqrt(3)r[2,1] -2x[2,2]             x[2,3]+sqrt(3)r[2,3]]
    M[3,:] = [x[3,1]+sqrt(3)r[3,1] x[3,2]-sqrt(3)r[3,2] -2x[3,3]            ]
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
            Pij = zeros(3,3)
            Qij = zeros(3,3)

            tr_vᵢ = value(variable_by_name(m, string("vⱼ", "[", i, "]")))
            Pij = value(variable_by_name(m, string("Pij", "[", i, "-", j, "]")))
            Qij = value(variable_by_name(m, string("Qij", "[", i, "-", j, "]")))
            
            Sᵢⱼ = @. complex(Pij, Qij)
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


"""
    reg_busses(p::Inputs)

All of the regulated busses, i.e. the second bus in the regulated edges
"""
function reg_busses(p::Inputs)
    getindex.(keys(p.regulators), 2)
end


function turn_ratio(p::Inputs, b::String)
    if !(b in reg_busses(p))
        throw(@error "Bus $b is not a regulated bus")
    end
    for (edge_tuple, d) in p.regulators
        if edge_tuple[2] == b
            return d[:turn_ratio]
        end
    end
end


function has_vreg(p::Inputs, b::String)
    for (edge_tuple, d) in p.regulators
        if edge_tuple[2] == b  && :vreg in keys(d)
            return true
        end
    end
    return false
end


function vreg(p::Inputs, b::String)
    if !(b in reg_busses(p))
        throw(@error "Bus $b is not a regulated bus")
    end
    for (edge_tuple, d) in p.regulators
        if edge_tuple[2] == b  && :vreg in keys(d)
            return d[:vreg]
        end
    end
    false
end