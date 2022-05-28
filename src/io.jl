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
    function dss_parse_lines(fp::String)

Parse a openDSS line codes file, returning 
    - edges, an array of tuples with 2 string values each for the sending and receiving busses on each edge
    - linecodes, an array of string
    - linelengths, an array of float

TODO is there a standard order for the line values? This function assumes that the order is:
    New Line.650632    Phases=3 Bus1=RG60.1.2.3   Bus2=632.1.2.3  LineCode=mtx601 Length=2000 units=ft
"""
function dss_parse_lines(fp::String)
    edges = Tuple[]
    linecodes = String[]
    linelengths = Float64[]
    for line in eachline(fp)
        if startswith(line, "New Line.") && !(contains(line, "Switch"))
            N = length(line)
            name = chop(line, head=findfirst(".", line)[end], tail=N-findfirst("  ", line)[1]+1)
            bus1 = strip(chop(line, head=findfirst("Bus1=", line)[end], tail=N-findfirst("Bus2", line)[1]+1))
            bus2 = strip(chop(line, head=findfirst("Bus2=", line)[end], tail=N-findfirst("LineCode", line)[1]+1))
            LineCode = strip(chop(line, head=findfirst("LineCode=", line)[end], tail=N-findfirst("Length", line)[1]+1))
            if contains(line, "units")
                LineLength = strip(chop(line, head=findfirst("Length=", line)[end], tail=N-findfirst("units", line)[1]+1))
            else
                LineLength = strip(chop(line, head=findfirst("Length=", line)[end], tail=0))
            end
            b1 = chop(bus1, tail=length(bus1)-findfirst(".", bus1)[1]+1)
            b2 = chop(bus2, tail=length(bus2)-findfirst(".", bus2)[1]+1)
            push!(edges, (b1, b2))
            push!(linecodes, convert(String, LineCode))
            push!(linelengths, parse(Float64, LineLength))
        end
    end
    return edges, linecodes, linelengths
end


function dss_parse_line_codes(fp::String, linecodes::Array{String, 1})
    d = Dict(c => Dict{String, Any}() for c in linecodes)
    open(fp) do io
        while !eof(io)
            line = readline(io)
            N = length(line)
            if startswith(line, "New linecode")
                code = convert(String, chop(line, head=findfirst(".", line)[end], tail=N-findfirst("nphases", line)[1]+2))
                if code in linecodes
                    while !occursin("rmatrix", line) || startswith(line, "!")
                        line = readline(io)
                    end
                    d[code]["rmatrix"] = dss_parse_string_matrix(line)

                    while !occursin("xmatrix", line) || startswith(line, "!")
                        line = readline(io)
                    end
                    d[code]["xmatrix"] = dss_parse_string_matrix(line)
                end
            end
        end
    end

    return d
end


function dss_parse_string_matrix(line::String)
    N = length(line)
    str_array = split(chop(line, head=findfirst("[", line)[end], tail=N-findfirst("]", line)[1]+1))
    filter!(s -> !occursin("|", s), str_array)
    return parse(Float64, str_array[1])
end