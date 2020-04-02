# Getting Started
* Installation
* Building a Model JuMP Model
* Solving
* Mathematical Options


## Installation
MaximinOPF.jl is availalbe to use as a local package with the Julia package manager. 
* [Download the repository](https://github.com/Argonne-National-Laboratory/MaximinOPF.jl/archive/master.zip)
* Install and test pacakge
```
(Shell)> cd <proj_root>
(Shell)> julia
(Julia)> ]
(v1.3) pkg> activate .
(MaximinOPF) pkg> test
... testing ...
Testing MaximinOPF tests passed
(MaximinOPF) pkg>
```

## Using a Package

### Module Importing
MaximinOPF.jl contains MaximinOPF module to build a JuMP model. The module provides two key [APIs](../API/) to solve the [mathematical problems](../mathematicalconcept/). The package uses [PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl/tree/master/docs) to parse the standard power system input files and to define the mathematical form of the model. 

```julia
using MaximinOPF
using PowerModels
```

### Input Preparing

In this example, we will use IEEE 9-BUS power system as a sample input. The input files are available on the data folder. The first step is that preparing the inputs for the MaximinOPF module with the PowerModels.

```julia
pm_data = PowerModels.parse_file(./data/case9.m)
powerform = SOCWRConicPowerModel
```

The pm_data is the ParsedPowerSystemModel object that describes the input power system. The powerform indicates the mathematical form of the model.

### Model Building

The building process is simply operated by calling the [API functions](../API/) in MaximinOPF module. In this example, we will call the MaximinOPFModel API to build a model.

```julia
maxmin_model = MaximinOPF.MaximinOPFModel(pm_data, powerform)
```

The output of the function is a JuMP model that includes the power flow, feasibility, minmax, dualization problems in SOCWRConic form. 

### Solving

In this example, we will use [Mosek.jl](https://github.com/JuliaOpt/Mosek.jl) as a solver interface. The solving code follows the used interfaces.
```julia
using Mosek
using MosekTools

set_optimizer(maxmin_model,Mosek.Optimizer)  
result = @elapsed JuMP.optimize!(maxmin_model)
#Print Result   
status=JuMP.termination_status(maxmin_model)
println("Time taken to solve is: ", result, " with status ",status,".")
```


## Mathematical Options
The MaximinOPF.jl uses [PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl/tree/master/docs) to build power system variables and constraints in different form of the optimization problem.
* TODO: Description of available PowerModel Options
