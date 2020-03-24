
# CAPS.jl

## Overview
CAPS.jl is a Julia/JuMP package for modeling N-k Contingency Analysis in Power Systems (CAPS), which creates a maxmin model where an attacker seeks to maximize system infeasibility subject to 1) a given attack budget of K network assets (typically lines and transformers) and 2) the assumed optimal defense response of the power system operator. The power flow system is modeled according to various formulations based on relaxations or approxiations of the equations governing line/transformer power flow, which include the SDP and SOCP relaxations, and the DC or other linear approximation of the power flow equations.

## Getting Started

* Installation
* Building a model
* Optimization options

### Installation
The CAPS.jl is available to use as a local package with the Julia package manager. 
```
(Shell)> cd <proj_root>
(Shell)> julia
(Julia)> ]
(v1.3) pkg> activate .
(MaximinOPF) pkg> test
```
### Building a model

### Optimization options
