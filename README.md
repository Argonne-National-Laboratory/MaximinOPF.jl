# MaximinOPF.jl

## Overview
MaximinOPF.jl is a Julia/JuMP package for modeling N-k Contingency Analysis in Power Systems (CAPS), which creates a maxmin model where an attacker seeks to maximize system infeasibility subject to 1) a given attack budget of K network assets (typically lines and transformers) and 2) the assumed optimal defense response of the power system operator. The power flow system is modeled according to various formulations based on relaxations or approxiations of the equations governing line/transformer power flow, which include the SDP and SOCP relaxations, and the DC or other linear approximation of the power flow equations.

## Quick Start
The simplest way to install the MaximinOPF.jl is below:
```
(Shell)> cd <proj_root>
(Shell)> julia
(Julia)> ]
(v1.3) pkg> activate .
(MaximinOPF.jl) pkg> build
(MaximinOPF.jl) pkg> test
```
The detail information to use this package is available on the [Getting Started](./gettingstarted/) document. The document include the installation, CAPS model building, and optimization option information.

## Links
* Code Repository: [MaximinOPF.jl](https://github.com/Argonne-National-Laboratory/MaximinOPF.jl)
* User Guide: [Getting Started](./gettingstarted/)
* API Documentation: [API Doc](./API/)
* Mathematical Concept: [Mathematical Concept](./mathematicalconcept/)

## Team Members
* Kibaek Kim: Mathematical Researcher, Project Leader
* Brian Dandurand: Mathematical Researcher
* Sang-il Yim: Software Developer

## License
* Argonne License (probably?)
* PowerModels License?
