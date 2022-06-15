module LinDistFlow

using JuMP
using LinearAlgebra
import SparseArrays: sparse
import PowerModelsDistribution: parse_dss

export 
    Inputs,
    singlephase38linesInputs,
    dsstxt_to_sparse_array, 
    parse_dss,
    build_ldf!,
    add_variables,
    constrain_power_balance,
    constrain_substation_voltage,
    constrain_KVL,
    constrain_loads,
    constrain_bounds,
    get_bus_values, 
    get_edge_values, 
    # recover_voltage_current,  # TODO validate this method
    i_to_j, 
    j_to_k, 
    rij, 
    xij

include("types.jl")
include("io.jl")
include("inputs.jl")
include("utils.jl")
include("model.jl")

end # module
