
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


function get_ijphases(i::String, j::String, p::Inputs{MultiPhase})
    ij_idx = get_ij_idx(i, j, p)
    return p.phases[ij_idx]
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
    remove_bus!(j::String, p::Inputs{MultiPhase})

Remove bus `j` in the line i->j->k from the model by making an equivalent line from busses i->k
"""
function remove_bus!(j::String, p::Inputs{MultiPhase})
    # TODO mv this to CommonOPF (by making rij & xij inputs to this method?)
    # get all the old values
    i, k = i_to_j(j, p)[1], j_to_k(j, p)[1]
    if j in reg_busses(p)
        for (edge, regdict) in p.regulators
            if edge[2] == j
                @warn "Moving regulator bus $j to bus $k in remove_bus!"
                p.regulators[(i,k)] = regdict
            end
        end
    end
    ij_idx, jk_idx = get_ij_idx(i, j, p), get_ij_idx(j, k, p)
    ij_len, jk_len = p.linelengths[ij_idx], p.linelengths[jk_idx]
    ij_linecode, jk_linecode = get_ijlinecode(i,j,p), get_ijlinecode(j,k,p)
    r_ij, x_ij, r_jk, x_jk = rij(i,j,p)*p.Zbase, xij(i,j,p)*p.Zbase, rij(j,k,p)*p.Zbase, xij(j,k,p)*p.Zbase
    phases = p.phases[ij_idx]
    # make the new values
    r_ik = r_ij .+ r_jk
    x_ik = x_ij .+ x_jk
    ik_len = ij_len + jk_len
    ik_linecode = ik_key = i * "-" * k
    ik_amps = minimum([p.Isquared_up_bounds[ij_linecode], p.Isquared_up_bounds[jk_linecode]])
    # delete the old values
    delete_edge_ij!(i, j, p)
    delete_edge_ij!(j, k, p)
    delete_bus_j!(j, p)
    # add the new values
    push!(p.edges, (i, k))
    push!(p.linecodes, ik_linecode)
    push!(p.phases, phases)
    push!(p.linelengths, ik_len)
    push!(p.edge_keys, ik_key)
    p.Zdict[ik_linecode] = Dict(
        "nphases" => length(phases),
        "name" => ik_linecode,
        "rmatrix" => r_ik ./ ik_len,
        "xmatrix" => x_ik ./ ik_len,
    )
    p.Isquared_up_bounds[ik_linecode] = ik_amps
end


"""
    combine_parallel_lines!(p::Inputs{MultiPhase})

Combine any parallel single phase lines without loads on intermediate busses into one multiphase line
"""
function combine_parallel_lines!(p::Inputs{MultiPhase})
    # TODO mv this to CommonOPF (by making rij & xij inputs to this method?)
    g = make_graph(p.busses, p.edges)
    end_bs = busses_with_multiple_inneighbors(g)

    for b2 in end_bs
        ins = inneighbors(g, b2)
        start_bs = unique(
            next_bus_above_with_outdegree_more_than_one.(repeat([g], length(ins)), ins)
        )
        if length(start_bs) == 1 && typeof(start_bs[1]) == String
            # we have a start bus and end bus to merge lines (if none of the intermediate busses have loads)
            b1 = start_bs[1]
            paths = paths_between(g, b1, b2)
            check_paths(paths, p)
            # remove all the intermdiate busses s.t. we have two // lines from b1 to b2
            for path in paths
                for b in path
                    remove_bus!(b, p)
                end
            end
            # now we combine the two // lines into one
            ij_idxs = findall(t->(t[1]==b1 && t[2]==b2), p.edges)
            @assert length(ij_idxs) in [2,3] "Found more than three parallel lines between busses $b1 and $b2!"
            if length(ij_idxs) == 2
                # old values
                i1, i2 = ij_idxs
                len1, len2 = p.linelengths[i1], p.linelengths[i2]
                phases1, phases2 = p.phases[i1], p.phases[i2]
                linecode1, linecode2 = p.linecodes[i1], p.linecodes[i2]
                rmatrix1, rmatrix2 = p.Zdict[linecode1]["rmatrix"] .* len1, p.Zdict[linecode2]["rmatrix"] .* len2
                xmatrix1, xmatrix2 = p.Zdict[linecode1]["xmatrix"] .* len1, p.Zdict[linecode2]["xmatrix"] .* len2
                amps1, amps2 = p.Isquared_up_bounds[linecode1], p.Isquared_up_bounds[linecode2]
                # new values
                new_len = (len1 + len2) / 2
                new_linecode = "combined__" * linecode1
                new_rmatrix = Diagonal([rmatrix1[1], rmatrix2[1]]) ./ new_len
                new_xmatrix = Diagonal([xmatrix1[1], xmatrix2[1]]) ./ new_len
                # TODO add off diagonal r/x values from up/downstream lines?
                new_phases = [phases1[1], phases2[1]]

                # delete the old values
                delete_edge_index!(i1, p)
                ij_idxs = findall(t->(t[1]==b1 && t[2]==b2), p.edges)
                delete_edge_index!(ij_idxs[1], p)
                
                # add the new values
                push!(p.edges, (b1, b2))
                push!(p.linecodes, new_linecode)
                push!(p.phases, new_phases)
                push!(p.linelengths, new_len)
                push!(p.edge_keys, b1 * "-" * b2)
                p.phases_into_bus[b2] = new_phases
                p.Zdict[new_linecode] = Dict(
                    "nphases" => 2,
                    "name" => new_linecode,
                    "rmatrix" => new_rmatrix,
                    "xmatrix" => new_xmatrix,
                )
                p.Isquared_up_bounds[new_linecode] = (amps1 + amps2) / 2
            else  # 3 lines to combine
                # old values
                i1, i2, i3 = ij_idxs
                len1, len2, len3 = p.linelengths[i1], p.linelengths[i2], p.linelengths[i3]
                phases1, phases2, phases3 = p.phases[i1], p.phases[i2], p.phases[i3]
                linecode1, linecode2, linecode3 = p.linecodes[i1], p.linecodes[i2], p.linecodes[i3]
                rmatrix1, rmatrix2, rmatrix3 = p.Zdict[linecode1]["rmatrix"] .* len1, p.Zdict[linecode2]["rmatrix"] .* len2, p.Zdict[linecode3]["rmatrix"] .* len3
                xmatrix1, xmatrix2, xmatrix3 = p.Zdict[linecode1]["xmatrix"] .* len1, p.Zdict[linecode2]["xmatrix"] .* len2, p.Zdict[linecode3]["xmatrix"] .* len3
                amps1, amps2, amps3 = p.Isquared_up_bounds[linecode1], p.Isquared_up_bounds[linecode2], p.Isquared_up_bounds[linecode3]
                # new values
                new_len = (len1 + len2 + len3) / 3
                new_linecode = "combined__" * linecode1 *"__"* linecode2 *"__"* linecode3
                new_rmatrix = Diagonal([rmatrix1[1], rmatrix2[1], rmatrix3[1]]) ./ new_len
                new_xmatrix = Diagonal([xmatrix1[1], xmatrix2[1], xmatrix3[1]]) ./ new_len
                # TODO add off diagonal r/x values from up/downstream lines?
                new_phases = [phases1[1], phases2[1], phases3[1]]

                # delete the old values
                delete_edge_index!(i1, p)
                ij_idxs = findall(t->(t[1]==b1 && t[2]==b2), p.edges)
                delete_edge_index!(ij_idxs[1], p)
                ij_idxs = findall(t->(t[1]==b1 && t[2]==b2), p.edges)
                delete_edge_index!(ij_idxs[1], p)
                
                # add the new values
                push!(p.edges, (b1, b2))
                push!(p.linecodes, new_linecode)
                push!(p.phases, new_phases)
                push!(p.linelengths, new_len)
                push!(p.edge_keys, b1 * "-" * b2)
                p.phases_into_bus[b2] = new_phases
                p.Zdict[new_linecode] = Dict(
                    "nphases" => 3,
                    "name" => new_linecode,
                    "rmatrix" => new_rmatrix,
                    "xmatrix" => new_xmatrix,
                )
                p.Isquared_up_bounds[new_linecode] = (amps1 + amps2 + amps3) / 3
            end
            @info "Made new combined line between busses $b1 and $b2"
        end
    end
end

