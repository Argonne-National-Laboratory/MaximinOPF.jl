The proximal bundle methods implemented here solve problems of the following form:

min c^T x + f(x)
s.t. Ax <= b
     g(x) <= 0

where 

1) c^T x constitute the known linear terms of the objective function, 
2) f(x) is a convex function that may be nonlinear, nonsmooth, and with unknown a priori closed-form
3) Ax <= b are the a priori-known affine constraints,
4) g(x) <= 0 is a convex constraint defined by a possibly nonlinear, nonsmooth function g(x) that is not of a priori-known closed form.
