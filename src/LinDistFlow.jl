module LinDistFlow

using CommonOPF
using JuMP
using LinearAlgebra
import SparseArrays: sparse

export 
    Inputs,
    singlephase38linesInputs,
    dsstxt_to_sparse_array, 
    dss_files_to_dict,
    build_ldf!,
    add_variables,
    constrain_power_balance,
    constrain_substation_voltage,
    constrain_KVL,
    constrain_loads,
    constrain_line_amps,
    get_bus_values, 
    get_edge_values, 
    # recover_voltage_current,  # TODO validate this method
    i_to_j, 
    j_to_k, 
    rij, 
    xij,
    reg_busses,
    turn_ratio,
    vreg,
    has_vreg,
    make_graph,
    remove_bus!,
    combine_parallel_lines!

include("utils.jl")
include("model_single_phase.jl")
include("model_multi_phase.jl")

end # module
