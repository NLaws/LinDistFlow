"""
    dsstxt_to_sparse_array(fp::String, first_data_row::Int = 5)

convert a SystemY.txt file from OpenDSS to a julia matrix.
assumes that Y is symmetric.
"""
function dsstxt_to_sparse_array(fp::String, first_data_row::Int = 5)

    rows = Int[]
    cols = Int[]
    real = Float64[]
    imag = Float64[]

    for (i, line) in enumerate(eachline(fp))

        if i < first_data_row continue end
        line = replace(line, " "=>"")  # "[1,1]=50500+j-50500"
        N = length(line)
        if N == 0 continue end

        append!(rows, tryparse(Int64,
                chop(line, head=findfirst("[", line)[end], tail=N-findfirst(",", line)[end]+1)
        ))

        append!(cols, tryparse(Int64,
                chop(line, head=findfirst(",", line)[end], tail=N-findfirst("]", line)[end]+1)
        ))

        append!(real, tryparse(Float64,
                chop(line, head=findfirst("=", line)[end], tail=N-findfirst("+", line)[end]+1)
        ))

        append!(imag, tryparse(Float64,
                chop(line, head=findfirst("j", line)[end], tail=0)
        ))
    end
    return convert(Array{Complex, 2}, Symmetric(sparse(rows, cols, complex.(real, imag)), :L))
end


"""
    dss_dict_to_arrays(d::Dict)

Parse the dict from PowerModelsDistribution.parse_dss into values needed for LinDistFlow
"""
function dss_dict_to_arrays(d::Dict)
    # TODO allocate empty arrays with number of lines
    edges = Tuple[]
    phases = Vector[]
    linecodes = String[]
    linelengths = Float64[]

    for v in values(d["line"])
        if "switch" in keys(v) && v["switch"] == true
            continue
        end
        try
            if occursin(".", v["bus1"])  # have to account for .1.2 phases for example
                b1 = chop(v["bus1"], tail=length(v["bus1"])-findfirst('.', v["bus1"])+1)
                phs = sort!(collect(parse(Int,ph) for ph in split(v["bus1"][findfirst('.', v["bus1"])+1:end], ".")))
            else  # default to 3 phases
                b1 = v["bus1"]
                phs = 123
            end

            if occursin(".", v["bus2"])
                b2 = chop(v["bus2"], tail=length(v["bus2"])-findfirst('.', v["bus2"])+1)
            else
                b2 = v["bus2"]
            end

            push!(edges, (b1, b2))
            push!(linecodes, v["linecode"])
            push!(linelengths, v["length"])
            push!(phases, phs)
        catch
            @warn("Unable to parse line $(v) when processing OpenDSS model.")
        end
    end

    return edges, linecodes, linelengths, d["linecode"], phases
end


"""
    dss_loads(d::Dict)

Return the P,Q loads from the dict provided by parse_dss. Indexed on bus (String), phase (integer),
then time (integer)

- TODO other reactive power specifications ?
- TODO handle vectors ("yearly" loads). Does the parse_dss handle redirects to txt files for loads?
- TODO LoadXfmrs.dss in 8500 Node ntwk: all loads are put on phases 1 or 2, 
    but actual MV phase could be 1, 2, or 3 -> have to map load phases through transformers?
"""
function dss_loads(d::Dict)
    P, Q = Dict{String, Dict{Int, Array{Real}}}(), Dict{String, Dict{Int, Array{Real}}}()
    for v in values(d["load"])
        bus = chop(v["bus1"], tail=length(v["bus1"])-findfirst('.', v["bus1"])+1)
        phases = collect(parse(Int,ph) for ph in split(v["bus1"][findfirst('.', v["bus1"])+1:end], "."))
        if !(bus in keys(P))
            P[bus] = Dict{Int, Array{Real}}()
            Q[bus] = Dict{Int, Array{Real}}()
        end
        if v["phases"] == 1
            P[bus][phases[1]] = [v["kw"] * 1000]  
            if "kvar" in keys(v)
                Q[bus][phases[1]] = [v["kvar"] * 1000]
            elseif "pf" in keys(v)
                p = P[bus][phases[1]]
                Q[bus][phases[1]] = sqrt( (p/v["pf"])^2 - p^2 )
            else
                Q[bus][phases[1]] = 0.0
            end
        else  # split the load evenly across phases
            p = v["kw"] / length(phases) * 1000
            
            if "kvar" in keys(v)
                q = v["kvar"] / length(phases) * 1000
            elseif "pf" in keys(v)
                q = sqrt( (p/v["pf"])^2 - p^2 )
            else
                q = 0.0
            end
            for phs in phases
                P[bus][phs] = [p]
                Q[bus][phs] = [q]
            end
        end
    end
    return P, Q
end
