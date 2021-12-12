# Linear-Dynamic-Response-History-Analysis-for-MDOF-systems
Linear dynamic response history analysis for estimation of structural response to dynamic loads (earthquakes, vibrations, shocks, etc.)

Matlab code for the application of the linear dynamic response history analysis (linear DRHA) of multi-degree of freedom (MDOF) structures is presented. Two procedures to calculate the dynamic response of a MDOF system subject to dynamic loading are included:
(a) direct integration of equations of motion and 
(b) the modal superposition of dynamic responses of SDOF systems equivalent to the eigenmodes to which a MDOF system is decomposed. 

For the direct integration of equations of motion, the function LDRHA_DI_MDOF.m is used. See the following examples:
•	example_Industrial_Building_DI_NPTEL.m
•	example_Shear_Building_2_DI_NPTEL.m
•	example_Shear_Frame_4A_DI_Chopra.m
•	example_Shear_Frame_5_DI_Chopra.m
for more details. 

The function LDRHA_DI_MDOF.m needs the damping matrix of the structure as user input. In the case that the damping of the structure is defined in terms of the critical damping ratios of its various eigenmodes and not by a damping matrix, the function CDM can be used to generate the classical damping matrix of the structure. Then this matrix can be input to the function LDRHA_DI_MDOF.m. See the following example:
•	example_Damping_Chopra
for more details.

For the modal superposition procedure for the dynamic response history analysis the functions LDRHA_MS_MDOF.m and LDRHA_SDOF.m are used. The latter is called inside the for-mer. See the following examples:
•	example_Shear_Frame_5_MS_Chopra.m
•	example_Industrial_Building_MS_NPTEL.m
•	example_Shear_Building_2_MS_NPTEL.m
•	example_Shear_Frame_4A_MS_Chopra.m
for more details. 

All the above functions can be used for acceleration time histories of a constant time step size. If this is not the case, then the acceleration time history needs to be resampled by using the MATLAB program file function RESAMPLE.m. The user is encouraged to see the example
•	example_Resampling_Nonuniform_Time_History.m
in this last case. 

The dynamic response history analysis procedure with direct integration proceeds incremental-ly, by solving the MDOF system equations for each time step (LDRHA_DI_MDOF.m, LDRHA_SDOF.m). 

The modal superposition dynamic response history analysis procedure (LDRHA_MS_MDOF.m) utilizes the following steps:
1. Define the structural properties.
a. Determine the mass matrix m and stiffness matrix k.
b. Estimate the modal damping ratios ζn.
2. Determine the natural frequencies ωn and modes φn.
3. Compute the response in each mode n by the following steps:
a. Compute the dynamic response qn(t) of a SDOF system with natural frequency ωn and damping ratio ζn.
b. Compute the nodal displacements from φn and qn(t)
c. Compute the element forces associated with the nodal displacements from φn, k and qn(t)
4. Combine the contributions of all the n modes to determine the total response.

The present code is accompanied by 10 examples in which its application is presented. These examples are taken from various standard textbooks or other material. The results of the ex-amples are verified by the results of the application of the present code.
The author is open to any suggestions or recommendations that the users may have.
