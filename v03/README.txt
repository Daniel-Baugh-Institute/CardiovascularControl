READme

Michelle Gee
Closed-loop modeling of intrinsic cardiac nervous system contributions to respiratory sinus arrhythmia
10/6/22

This folder includes scripts to produce Figures 3-7 and the tables of the manuscript. Each script can be run on its own in MATLAB to produce each figure.

Specific figure notes

Figure 5: The data is provided in RSA_GPR_kRSA_09_30_22.mat so Fig5_kRSA_likelihood can be run to produce Figure 5B. To reproduce the data in the .mat file, files in the folder Cluster_required_Fig5_paramTuning were used. After connecting to an HPC cluster using a slurm job handler, the batch script pbatch_Fig5_GPR.qs was submitted using the command sbatch pbatch_Fig5_GPR.qs. This runs the script Fig5_GPR_cluster.m and the associated Simulink model ICN_model_v15_VNS_cluster.slx. The directory Cluster_required_Fig5_paramTuning is designed such that the entire directory should be transferred to the cluster so that the relevant functions are available.

Figure 7: To download the experimental data, go to https://sparc.science/datasets/28?type=dataset&datasetDetailsTab=files&path=files%2Fprimary and download Pig013_ICNS15_Matlab.mat. Add it to the v03 directory with the rest of the files. The figure can be produced by running Fig7_VNS.m because the data are provided in Aff1_Eff0_5.mat, AfferentOnly.mat, and Both1Hz.mat. The model is set to produce results for a 1 Hz afferent stimulation and 0.5 Hz efferent stimulation. The stimulation parameters can be adjusted in the simulink model in the Autonomic Nervous System block by adjusting the magnitude of the blocks Step4-7 to match the desired stimulation parameters.

ICN parameter tuning: ICN parameter tuning can be done by submitting pbatch_ICNtune.qs to a cluster, which then runs the script DARWIN_BatchRunParallel_1558902.m and the associated Simulink files ICN_model_v15_tune.slx. The scrpit is currently set to reproduce the ICN parameter values selected for this manuscript. The output will be in the form of a .out file, which can be opened with any text editor software. The file contains the top 10 parameter sets (in order from best to worst) from two rounds of Sobol sampling selected for the heart rate, elastance, and ICN neuron firing frequency metrics that lead to the lowest mean square error compared to experimental values (Rajendran 2019). Firing frequency values from Rajendran 2019 were calculated using the provided script Analysis_Rajendran.m.

Rajendran, P., Vaseghi, M., & Ardell, J. (2019). Functional recordings from the pig intrinsic cardiac nervous system (ICN) (Version 2) [Data set]. SPARC Consortium. https://doi.org/10.26275/OWRI-MPSX
