Modeling and Analysis of the Intrinsic Cardiac Nervous System in Cardiovascular Control and Respiratory Sinus Arrhythmia

Michelle Gee
July 18 2022

Figures 2-5 were created by running files Fig2_RSA_Validation.m, Fig3_4_HRandEmax_perturb_Suga.m, and Fig5_plot.m on Matlab R2020b (Simulink will also be necessary). The associated Simulink simulation files are included in the folder and will run when the associated script is executed.

Figures 6 and 7 were generated using the University of Delaware's DARWIN computing resources to run simulations in parallel. Files for these figures are included in the subfolder Fig6_7_Cluster_Required. A submission script pbatch_Fig6.qs (Figure 6) or pbatch_PAWN.qs (Figure 7), was submitted to the job scheduler (slurm) used by University of Delaware. The submission script calls a Matlab script RSA_GPR_robuts_v3.m (Figure 6) or PAWN_CardioControl_v3.m (Figure 7) that runs the simulations in parallel and saves the output files in .mat format. The results can then be plotted and analyzed using Fig6_plot.m or Fig7_plot.m. If parallel computing resources are not available, the .mat files have been provided and the plotting scripts Fig6_plot.m and Fig7_plot.m can be run independently.