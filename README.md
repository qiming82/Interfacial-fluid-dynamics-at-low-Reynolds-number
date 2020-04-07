# Boundary element method 

Boundary integral code to compute a perfectly conducting core-annular (axisymmetric case two phase flow at low Reynolds number) viscous thread subject to a radial electric field.

The code only works for viscosity ratio \lambda = 1 but qualitatively captured the behavior of the system when \lambda is close to one. Code attached here is a test version and no guarantee for your own problem.

One of the difficult parts is the computation of the periodic (axisymmetric) Green's functions for both Laplace and Stokes flow problem. The periodic Green's function for Stokes flow problem is adapted from the code by C. Pozrikidis (along with some simple auxiliary functions from Pozrikidis http://dehesa.sourceforge.net/BEMLIB/), whereas for Laplace equation, the implementation is similar (the periodic version) as the leading order singulairy of Green's function for Laplace equation is the same as the Stokeslet (both weakly singular).

Results to see or cite the Ref. paper if using the code: 
"Qiming Wang and Demetrios T. Papageorgiou, Dynamics of a viscous thread surrounded by another viscous fluid in a cylindrical tube under the action of a radial electric field: Breakup and touchdown singularities J. Fluid Mech. 2011, Vol 683, 27-56."
