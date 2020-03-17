# CAPS.jl APIs
CAPS.jl provides two key APIs to build a JuMP model based on the [PowerModels.jl's](https://github.com/lanl-ansi/PowerModels.jl/tree/master/docs) ParsedPowerSystemModel.

## MaximinOPFModel(ParsedPowerSystemModel, MathOption)
This API build a JuMP model for the [maximin]() problem that contains [minmax]() and [dualization]() problems.
The code examples:

```julia
using MaximinOPF
using PowerModels

pm_data = PowerModels.parse_file(./data/case9.m)
powerform = SOCWRConicPowerModel

maxmin_model = MaximinOPF.MaximinOPFModel(pm_data, powerform)
```
	
The sample code uses IEEE 9-bus test system as a input power system. PowerModels parses the system and builds a ParsedPowerSystemModel object. The SOCWRConicPowerModel, which is defined in the PowerModels, option indicates the mathematical form of the model. The MaximinOPFModel function build a JuMP model that can be solved by mathmatical solvers.

## MinimaxOPFModel(ParsedPowerSystemModel, MathOption)
This API build a power system model for the [minmax problem]() that contains [feasibility]() problem
The code examples:

```julia
using MaximinOPF
using PowerModels

pm_data = PowerModels.parse_file(./data/case9.m)
powerform = SOCWRConicPowerModel

maxmin_model = MaximinOPF.MinimaxOPFModel(pm_data, powerform)
```

The sample code uses IEEE 9-bus test system as a input power system. PowerModels parses the system and builds a ParsedPowerSystemModel object. The SOCWRConicPowerModel, which is defined in the PowerModels, option indicates the mathematical form of the model. The MinimaxOPFModel function build a JuMP model that can be solved by mathmatical solvers.