# LinDistFlow Changelog

## v0.5.3
- upgrade CommonOPF to v0.3.8

## v0.5.2
- upgrade CommonOPF to v0.3.7

## v0.5.1
- add `get_line_amps`, `get_peak_line_amps_percent`, and `define_line_amps_pu`
- mv SinglePhase `rij` and `xij` to CommonOPF

## v0.5.0
- deprecate `get_edge_values` and `get_bus_values` in favor of `CommonOPF.get_variable_values`
- add `Results` for `SinglePhase` and `MultiPhase` models

## v0.4.1
- store `m[:injectioncons]` and `m[:loadbalcons]` in single phase models

## v0.4.0
- add single phase voltage regulators
- multiphase regulators index `:turn_ratio` on phase (breaking from CommonOPF)
- add `combine_parallel_lines!` and `remove_bus!`

## v0.3.2
- deprecate `constrain_bounds`
- mv `reg_busses`, `turn_ratio`, `has_vreg`, `vreg` to CommonOPF
-  split model.jl into model_single_phase.jl and model_multi_phase.jl

## v0.3.1
- update to CommonOPF v0.2.1

## v0.3.0
- move io.jl, types.jl, inputs.jl and some utils.jl to CommonOPF (new dependency)

## v0.2.0
- add three phase, unbalanced capability and documentation
- remove `dss_parse_lines` and `dss_parse_line_codes` methods; add `PowerModelsDistribution.parse_dss`
- `Inputs` only requires a "master.dss" file path now (instead two file paths for line codes and lines)

## v0.1.4
- update dependencies (JuMP 1.0, HiGHS for testing)

## v0.1.3
- remove CSV dependency

## v0.1.2
- upgrade Manifest.toml and Project.toml for Julia 1.7

## v0.1.1
- add `singlephase38linesInputs` type for quickly building a demo network
- deprecate `constraint_bounds` in favor of passing bounds to `Inputs` and constraining in `@variable`
