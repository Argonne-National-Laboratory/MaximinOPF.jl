
# CAPS.jl

## Overview
CAPS.jl is a Julia/JuMP package for modeling N-k Contingency Analysis in Power Systems. This package creates a maximin optimization model that maximizes the minimum system infeasibility by modeling an attacker and a defender of the power system operator. For a given budget K, the attacker finds a set of critical assets that maximizes the total system infeasibility to the system, whereas the defender minimizes the damange from the attack. The formulations supported for the underlying power system include the SDP and SOCP relaxations and DC approximation of the optimal power flow.

The proximal bundle methods implemented here solve problems of the following form:

$$
min c^T x + f(x) s.t. Ax <= b g(x) <= 0

where

c^T x constitute the known linear terms of the objective function,
f(x) is a convex function that may be nonlinear, nonsmooth, and with unknown a priori closed-form
Ax <= b are the a priori-known affine constraints
g(x) <= 0 is a convex constraint defined by a possibly nonlinear, nonsmooth function g(x) that is not of a priori-known closed form.
We assume for the sake of nontriviality that f and g are nonsmooth, and cannot be represented as the maximum value of a finite number of differentiable functions. The solution approach exchanges the solution of the above problem with iterations k of the following proximal point approximation centered about x^k

min c^T x + F^k(x) + (1/2t^k) * ||x-x^k||^2 s.t. Ax <= b G^k(x) <= 0 where F and G are cutting plane models of f and g, respectively.
$$
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
(v1.2) pkg> activate .
(MaximinOPF) pkg> test
```
### Building a model

### Optimization options
