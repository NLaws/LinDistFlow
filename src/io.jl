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
    # TODO ignore swithces?
    edges = Tuple[]
    linecodes = String[]
    linelengths = Float64[]

    for v in values(d["line"])
        b1 = chop(v["bus1"], tail=length(v["bus1"])-findfirst(".", v["bus1"])[1]+1)
        b2 = chop(v["bus2"], tail=length(v["bus2"])-findfirst(".", v["bus2"])[1]+1)
        push!(edges, (b1, b2))
        push!(linecodes, v["linecode"])
        push!(linelengths, v["length"])
    end

    return edges, linecodes, linelengths, d["linecode"]
end