SBEs solver 1.0
This is a script to solve simple, 1dimension, 2bands semiconduction bloch equations omiting Coulomb interaction, dephasing of occupation, background polarization, Berry phase.

===========================================
How to use this script?
a. Prepare the file "transition dipole moment.txt", which is corresponding to transition dipole moment of k grids.
b. Modify files "input_parameters.m" according to your input parameters.
c. Run file "sbe_solver.m"

ps:
a. After runnung this script successfully, you can get variables named "sol", which is saved in file named "sol.mat". Because the forth part of script is very time-consuming, so you can avoid to running this part again if you had run it  and gotten "sol.mat".
b. We fix the carriers occupation at boundary of Brillouin zone, it could cause error, so we suggest you take more k points(more than 50).
c. To improved: because it's a 1D model, youll overestimate the excited carriers along the laser field compared with some references of 3D model. Thus, we suggest you consider the carriers occupation(figure7 after running) relative value.
d. We use FFT of the current representing the spectrum rather than the first derivative of current with respect to time.  
e. Welcome to contact with me, email: vanking20000918@gmail.com.

============================================
Reference
[1] KIRA M, KOCH S W. Semiconductor quantum optics[M]. Cambridge University Press, 2011.
[2] GOLDE D, MEIER T, KOCH S W. High harmonics generated in semiconductor nanostructures by the coupled dynamics of optical inter- and intraband excitations[J/OL]. Physical Review B, 2008, 77(7). DOI:10.1103/physrevb.77.075330.
[3] GARG M, ZHAN M, LUU T T, et al. Multi-petahertz electronic metrology[J/OL]. Nature, 2016, 538(7625): 359-363. DOI:10.1038/nature19821.

• Copyright ©2023 VanKing, all rights reserved.

