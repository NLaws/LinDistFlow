# LinDistFlow Changelog

## dev
- deprecate `constrain_bounds`

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
