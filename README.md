
# CAPS.jl

## Overview
CAPS.jl is a Julia/JuMP package for modeling N-k Contingency Analysis in Power Systems. This package creates a maximin optimization model that maximizes the minimum system infeasibility by modeling an attacker and a defender of the power system operator. For a given budget K, the attacker finds a set of critical assets that maximizes the total system infeasibility to the system, whereas the defender minimizes the damage from the attack. The formulations supported for the underlying power system include the SDP and SOCP relaxations and DC approximation of the optimal power flow.

## Getting Started

* Installation
* Building a model
* Optimization options

### Installation
The CAPS.jl is availalbe to use as a local package with the Julia package manager. 
```
(Shell)> cd <proj_root>
(Shell)> julia
(Julia)> ]
(v1.3) pkg> activate .
(MaximinOPF) pkg> test
```
### Building a model

### Optimization options
