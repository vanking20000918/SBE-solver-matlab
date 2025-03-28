SBEs solver 1.0
This is a script to solve simple, 1dimension, 2bands semiconduction bloch equations omiting Coulomb interaction, dephasing of occupation, background polarization, Berry phase.
1D-2bands-semiconductor bloch equations:
iℏ ∂/∂t p_k=(ε_k^e+ε_k^h-i ℏ/T_2 ) p_k-(1-f_k^e-f_k^h ) Ω_k (t)+ieE(t)⋅∇_k p_k
ℏ ∂/∂t f_k^e(h) =-2 Im[Ω_k (t) p_k^* ]+eE(t)∙∇_k f_k^e(h) 
where p_k is micropolarization of k wave vector，ε_k^e is single particle energy of k wave vector on conduction band(e)，ε_k^h is single particle energy of k wave vector on valence band(h)，T2 is dephasing time，f_k^(e(h)) is carrier occupation on e(h)，Ω_k (t)=d_k E(t) is rabi frequency，and d_k is transition dipole moment，E(t) is electric field，ℏ,e is reduced planck constant amd eletmentary charge, respectively.
P(t)=∑_k[d_k p_k (t)+c.c.] 
P(t) is macroplolarization。
Interband current: J_er (t)=dP(t)/dt
Intraband current: J_ra (t)=e/ℏ ∑_k〖(∇_k E_k^e ) f_k^e 〗-e/ℏ ∑_k〖(∇_k E_k^h ) f_k^h 〗
Total current: J_tot (t)=J_ra (t)+J_er (t)
The above "E_k^(e(h))" expresses the band energy of e(h), which equals to single particle energy：ε_k^e=E_k^e,ε_k^e=〖-1⋅E〗_k^h.
Eventually, you could obtain the HHG with fouriel transform:
(I_HHG (ω)∝|FT{J_tot (t)} |^2
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

