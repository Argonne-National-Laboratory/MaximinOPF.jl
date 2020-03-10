# CAPS.jl

## Overview
CAPS.jl is a Julia/JuMP package for modeling N-k Contingency Analysis in Power Systems. This package creates a maximin optimization model that maximizes the minimum system infeasibility by modeling an attacker and a defender of the power system operator. For a given budget K, the attacker finds a set of critical assets that maximizes the total system infeasibility to the system, whereas the defender minimizes the damange from the attack. The formulations supported for the underlying power system include the SDP and SOCP relaxations and DC approximation of the optimal power flow.

The proximal bundle methods implemented here solve problems of the following form: 

$min(\mathcal{C}^\mathcal{T}\mathcal{x} + \mathcal{f}(x)) s.t. \mathcal{A}_\mathcal{x} <= \mathcal{b}, \mathcal{g}(x) <= 0$

where

* The $\mathcal{C}^\mathcal{T}\mathcal{x}$constitute the known linear terms of the objective function
* The $\mathcal{f}(x)$is a convex function that may be nonlinear, nonsmooth, and with unknown a priori closed-form
* The $\mathcal{A}_\mathcal{x} <= \mathcal{b}$are the a priori-known affine constraints
* The $\mathcal{g}(x) <= 0$is a convex constraint defined by a possibly nonlinear, nonsmooth function g(x) that is not of a priori-known closed form.

We assume for the sake of nontriviality that f and g are nonsmooth, and cannot be represented as the maximum value of a finite number of differentiable functions. The solution approach exchanges the solution of the above problem with iterations $\mathcal{k}$ of the following proximal point approximation centered about $\mathcal{x}^\mathcal{k}$

$min(\mathcal{C}^\mathcal{T}\mathcal{x} + \mathcal{F}^{\mathcal{k}(x)} + 1/2\mathcal{t}^\mathcal{k} * \mathopen||\mathcal{x}-\mathcal{x}^\mathcal{k}\mathclose||^2) s.t. \mathcal{A}_\mathcal{x} <= \mathcal{b}, $\mathcal{G}^{\mathcal{k}(x)} <= 0$

where
* The $\mathcal{F}$and $\mathcal{G}$ are cutting plane models of $\mathcal{f}$ and $\mathcal{g}$, respectively.

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
