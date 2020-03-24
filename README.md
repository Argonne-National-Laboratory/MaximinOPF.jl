# CAPS.jl

## Overview
CAPS.jl is a Julia/JuMP package for modeling N-k Contingency Analysis in Power Systems. This package creates a maximin optimization model that maximizes the minimum system infeasibility by modeling an attacker and a defender of the power system operator. For a given budget K, the attacker finds a set of critical assets that maximizes the total system infeasibility to the system, whereas the defender minimizes the damange from the attack. The formulations supported for the underlying power system include the SDP and SOCP relaxations and DC approximation of the optimal power flow.

## Quick Start
The simplest way to install the CAPS.jl is below:
```
(Shell)> cd <proj_root>
(Shell)> julia
(Julia)> ]
(v1.2) pkg> activate .
(CAPS.jl) pkg> build
(CAPS.jl) pkg> test
```
The detail information to use this package is available on the [Getting Started](./gettingstarted/) document. The document include the installation, CAPS model building, and optimization option information.

## Links
* Code Repository: [CAPS.jl](https://github.com/kibaekkim/CAPS.jl)
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
