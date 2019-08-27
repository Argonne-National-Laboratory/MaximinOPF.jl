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
