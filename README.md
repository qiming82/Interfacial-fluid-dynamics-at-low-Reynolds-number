# BIE_ecaf

Boundary integral code to compute a perfectly conducting core-annular viscous thread subject to a radial electric field.

The code only works for viscosity ratio \lambda = 1 but qualitatively captured the behavior of the system.

One of the difficult parts is the computation of the periodic (axisymmetric) Green's functions for both Laplace and Stokes flow problem. The periodic Green's function for Stokes flow problem is adapted from the code by C. Pozrikidis, whereas for Laplace equation, the implementation is similar.

Reference paper: 
"Qiming Wang and Demetrios T. Papageorgiou, Dynamics of a viscous thread surrounded by another viscous fluid in a cylindrical tube under the action of a radial electric field: Breakup and touchdown singularities J. Fluid Mech. 2011, Vol 683, 27-56."
