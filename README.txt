Multiple Fluid Prototype Solver 

Original SPH implementation implemented by Dongli Zhang (dongli.zhang0129@gmail.com) which can be found here: https://github.com/finallyjustice/sphfluid 

Modified by Edina Morris for multiple fluid simulation based on the paper "Multiple-fluid SPH Simulation Using a Mixture Model" by Ren et al. (Available at http://gamma.cs.unc.edu/SPH_MULTI_FLUIDS/TOG2014.pdf)

Controls:

1 = default scenario - side by side
2 = damn scenario
3 = drop scenario

Hit space bar to start simulation

Use arrows and WSAD to move camera around world and zoom, use mouse to change camera rotation

To change number of phases(liquids) you must alter the define in SphData.h and add individual phase data for the new liquid inside SphPhaseData.cpp