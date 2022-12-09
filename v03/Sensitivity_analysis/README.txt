PAWN sensitivity analysis using methodology and package from Pianosi and Wagener 2018.
The original package can be found here: https://www.safetoolbox.info/pawn-method/
The functions from the PAWN toolbox were copied over for the usef's convenience and are not presented as our original work. 

To run the sensitivity analysis, University of Delaware's DARWIN HPC was used to run the simulations. Skip to step 3 if a cluster is not available and you just want to plot the results of the sensitivity analysis.
1. Upload the entire sensitivity analysis directory to the HPC
2. Submit pbatch_PAWN.qs to the cluster using the sbatch command. This will run the script PAWN_analysis_main.m.
3. The output of PAWN_analysis_main.m includes PAWN_HR_v15.mat, which is the result of the sensitivity analysis. These results can be visualized by running Supplementary_Fig1.m on a local machine once PAWN_HR_v15.mat is in the same directory. 
