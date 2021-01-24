# Interfacial-fluid-dynamics-at-low-Reynolds-number
1. in folder 'drop', the code is used to 
    * simulate capillary pinchoff of viscous droplet immersed in another viscous fluid
    * use adaptive regridding to approach self-similar singularities before pinching (only for moderate viscosity ratio)
    * in principle, the code is able to reproduce the results in the follow two papers; if you find it helpful, please kindly cite the follow papers.
    * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

@article{BPSW2013,
Author = {M. R. Booty, D. T. Papageorgiou, M. Siegel and Q. Wang},
Title = {Long-wave equations and direct numerical simulations for the breakup of a viscous fluid thread surrounded by an immiscible viscous fluid},
Journal  = {IMA J. Applied Maths.},
volume = {78},
pages = {851-867},
Year = {2013}
}

@article{FW2021,
Author = {M. Fontelos and Q. Wang},
Title = {Discrete self-similarity in the formation of satellites for viscous cavity break up},
Journal  = {Phys. Rev. Fluids},
volume = {6},
pages = {013201}
Year = {2021}
}



2. in folder 'caff', the code is used to simulate conducting core annular flow in Stokes flow
   * the code handles computation of periodic Stokeslets and Laplace single layer Green's functions in axisymmetric geometry
   * the physical problem involves a core annular flow with equal viscosity where the core fluid can be a perfectly conducting liquid
   * if you find it helpful, please kindly cite the following paper
   * again, This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
   
@article{WP2011,
Author = {Q. Wang and D. T. Papageorgiou},
Title = {Dynamics of a viscous thread surrounded by another viscous fluid in a cylindrical tube under the action of a radial electric field: Breakup and touchdown singularities},
Journal  = {J. Fluid Mech.},
volume = {683},
pages = {27-56}
Year = {2011}
}
