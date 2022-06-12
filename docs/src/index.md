# LinDistFlow.jl

`LinDistFlow` builds the linear distflow constraints using JuMP. Currently only the single phase model is supported, but the multiphase unbalanced model is under progress.

The intent of this package is to allow users to build mathematical programs that include LinDistFlow constraints.
No objective is added to the JuMP model in this package and so solving any problem defined by the constraints built by LinDistFlow.jl is a feasibility problem. Dictionaries of constraints are provided so that one can delete and/or modify the base constraints to fit their problem.

# Inputs
There are two methods for creating `Inputs`:
1. Using openDSS files
2. Providing the network topology
```@docs
Inputs(::String, ::String)
Inputs(::Array{Tuple}, ::Array{String}, ::Array{Float64}, ::Vector{Vector}, ::String)
```
Both of the `Inputs` functions return a mutable `Inputs` struct:
```@docs
Inputs
```

# Building a Model
The `build_ldf!` function takes a `JuMP.Model` and `Inputs` struct as its two arguments and adds the variables and constraints:
```@docs
build_ldf!
```

# Variables
Let `m` be the JuMP.Model provided by the user, then the variables can be accessed via:
- `m[:vsqrd]` voltage magnitude squared, indexed on busses, (phases), time
- `m[:Pⱼ], m[:Qⱼ]` net real, reactive power injection, indexed on busses, (phases), time
- `m[:Pᵢⱼ], m[:Qᵢⱼ]` net real, reactive line flow, indexed on edges, (phases), time
After a model has been solved using `JuMP.optimize!` variable values can be extracted with `JuMP.value`. For more see [Getting started with JuMP](https://jump.dev/JuMP.jl/stable/tutorials/getting_started/getting_started_with_JuMP/#Getting-started-with-JuMP).


!!! note
    Single phase models do not have a phase index

# Accessing and Modifying Constraints
Let the JuMP.Model provided by the user be called `m`. All constraints are stored in `m[:cons]` as anonymous constraints.

## Power Injections
LinDistFlow.jl uses the convention that power injections are positive (and loads are negative). If no load is provided for a given bus (and phase) then the real and reactive power injections at that bus (and phase) are set to zero with an equality constraint.

All power injection constraints are stored in `m[:cons][:injection_equalities]`. The constraints are indexed in the following order:
1. by bus name (string), as provided in `Inputs.busses`;
2. by `:P` or `:Q` for real and reactive power respectively;
3. by phase number (integer); and
4. by time (integer).
For example, `m[:cons][:injection_equalities]["680"][:P][2][1]` contains the constraint reference for the power injection equality constraint for bus "680", real power, on phase 2, in time step 1.

If one wished to replace any constraint one must first delete the constraint using the `delete` function. For example:
```julia
delete(m, m[:cons][:injection_equalities]["680"][:P][1])
```
Note that the time index was not provided in the `delete` command in this example, which implies that the equality constraints for all time steps were deleted. One can also delete individual time step constraints by providing the time index.

The deleted constraints can then be replaced with a new set of constraints. For example:
```julia
m[:cons][:injection_equalities]["680"][:P][1] = @constraint(m, [t in 1:p.Ntimesteps],
    m[:Pⱼ]["680",1,t] == -1e3 / p.Sbase
)
```
where `p` is short for "parameters" and is the `Inputs` struct for the problem of interest. Note that it is not necessary to store the new constraints in the `m[:cons][:injection_equalities]`.

See the [JuMP documentation](https://jump.dev/JuMP.jl/stable/manual/constraints/#Delete-a-constraint) for more on deleting constraints.