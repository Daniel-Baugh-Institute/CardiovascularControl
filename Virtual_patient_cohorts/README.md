# Code reproduces simulations and figures from "Baroreflex adaptation following myocardial infarction in an in silico patient cohort"
# Step by step instructions are included for reproducing the figures below. Note that many of the simulations were run on a high performance computer (HPC) so the data files are provided so that the figures can be reproduced. The code to run the original simulations is included, but commented out.

# If you have any questions or find bugs please reach out to mmgee@udel.edu

# Reproducing the figures

# Fig 3A: in command line run plot_PAWN_scatter('PAWN_MAP_20240510.mat','PAWN_HR_20240510.mat','sensitivity_analysis.png'). Note that the PAWN toolbox (https://safetoolbox.github.io/Pawn.html) must be downloaded to run the full sensitivity analysis code using PAWN_MI_051024.m. 
# Fig 3B and C: Fig3_preIRmodels_060724.m
# Fig 3C: preIR_boxplots.m

# Fig 4C, D, E: postIR_main.m

# Fig 5B, C, D: postIR_baroreceptors.m
# Fig 5E: mdlParamstSNE_062324.m

# Fig 6B: overlay plots from postIR_NADMV_main.m and postIR_NTS_main.m
# Fig 6C and D: postIR_NTS_main.m
# Fig 6E and F: postIR_NADMV_main.m

# Fig 7B, C, D: postIR_ICN_main.m
# Fig 7E: central_peripheral_contributions.m

# Fig 8B, C, D: postIR_all_main.m

# Fig 9A, B, C: mdlParamstSNE_062324.m

# Producing and filtering the preIR population
# Run barocurve_dist_main.m
# This was run on an HPC (DARWIN HPC at University of Delaware). An example batch submission script is included (Mastitskaya_batch.qs).

# Producing and filtering the postIR populations
# postIR_XX_main.m, where XX could be baroreceptors, NTS, NADMV, ICN, all, or empty (cardiac). Again this was run on an HPC (DARWIN HPC at University of Delaware). An example batch submission script is included (Mastitskaya_batch.qs).

