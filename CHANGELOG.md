# LinDistFlow Changelog

## v0.1.4
- update dependencies (JuMP 1.0, HiGHS for testing)

## v0.1.3
- remove CSV dependency

## v0.1.2
- upgrade Manifest.toml and Project.toml for Julia 1.7

## v0.1.1
- add `singlephase38linesInputs` type for quickly building a demo network
- deprecate `constraint_bounds` in favor of passing bounds to `Inputs` and constraining in `@variable`
