The proximal bundle methods implemented here solve problems of the following form:

min c^T x + f(x)
s.t. Ax <= b
     g(x) <= 0

where 

<ol>
<li> c^T x constitute the known linear terms of the objective function, </li>
<li> f(x) is a convex function that may be nonlinear, nonsmooth, and with unknown a priori closed-form </li>
<li> Ax <= b are the a priori-known affine constraints </li>
<li> g(x) <= 0 is a convex constraint defined by a possibly nonlinear, nonsmooth function g(x) that is not of a priori-known closed form. </li>
</ol>

We assume for the sake of nontriviality that f and g are nonsmooth, and cannot be represented as the maximum value of a finite number of differentiable functions.
The solution approach exchanges the solution of the above problem with iterations k of the following proximal point approximation centered about x^k

min c^T x + F^k(x) + (1/2t^k) * ||x-x^k||^2
s.t. Ax <= b
     G^k(x) <= 0
where F and G are cutting plane models of f and g, respectively.


...

To run (on moonshot)

> ./runMoonshotProxPtSDP CASE K H

where CASE is the case instance, K is a nonnegative integer budget, and H=0,1,2 encodes a heuristic, with H=0 indicating no heuristic.
For example,

> ./runMoonshotProxPtSDP 30 4 0

runs the Case 30 instance with an attacker budget of 4 lines, and no heuristic. For more details on what is being done, see Main.jl and ProxPtSDP.jl
