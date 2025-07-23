# Mu-LaB-SWE
A multilayer shallow water equations solver using the Lattice Boltzmann method as part of a Masters research project. The module is coded in Fortran 90 and heavily inspired J.G. Zhou's works as well as Z. Ru et al. 

The module is still under development. For the moment, it can only simulate single layer simulations that have already been demonstrated effectively. The ultimate objective is to have a  multilayer shallow water equations and sediment transport solver using the Lattice Boltzmann Method with a turbulence model.

## Why the Lattice Boltzmann Method?
The Lattice Boltzmann Method (LBM) is a relatively recent method in computational hydrodynamics. However, it has been shown to be effective in subcritical (Fr < 1) regimes. Most applications of sediment transport are subcritical, with the exception of swash zones. What makes this method effective is its easy to parallelize. In an era of constantly upgrading GPU's, the methods which can take best advantage of parellization will be the most effective. It can also simulate around complex boundary conditions accurately, such as internal boundaries. Furthermore, it has a low diffusivity and simple, explicit computation. This high computing speed is compensated by its higher memory usage than other usual methods, like the Finite Volume Method (FVM). 
