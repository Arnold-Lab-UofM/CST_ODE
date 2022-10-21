## CST_ODE: MATLAB Scripts for "Evaluation of vaginal microbiome equilibrium states identifies microbial parameters linked to resilience after menses and antibiotic therapy"

Instructions to run code:
- Install Symbolic Math Toolbox, Parallel Computing Toolbox and Bioinformatics Toolbox
- Add ODE_scripts/ and Global_Sensitivity/ folders to your MATLAB Path
- Use Run_Files/ folders and MATLAB scripts to guide you through how to generate plots for eaach figure
    - workspaces/ folder provides workspaces generated and used for downstream analysis of our work

### Folder Layout

#### Run Files
The Run Files folder has scripts for generating each figure in the main text of the manuscript and the workspaces with any additional data needed to complete the analysis.

- Figure1
    - Fig1_run_global_SA.m
    - lhs_settings_input.xlsx _(input file for LHS pipeline)_
- Figure2
   - Fig2_compare_1SS_2SS.m
- Figure3
    - Fig3_run_Global_2D_Bifurcation_Plots.m
    - Fig3_run_LHS_analysis_menses.m
- Figure4
    - Fig4_run_1D_Global_Bifurcation.m
    - Fig4_run_LHS_analysis_abx.m
- Figure5
    - Fig5_run_LHS_analysis_abx_prebiotic.m
    - Fig5_run_LHS_analysis_BV_abx_dose_dur.m
- Workspaces
    - 2022-04-29-VALENCIA-16s-to-SS-type.mat _(Patient data by VALENCIA CSTs)_
    - HMP-UAB-MetaData.mat _(patient data on menses, antibiotic use, etc.)_
    - HMP-UAB-RelAbundance-Trajectories.mat _(clinical data - relative abundance over time)_
    - SSConfig-Analysis-Model_LHS10x.mat _(Equilvalent output from Figure1/ run script)_

#### ODE Scripts
The ODE Scripts folder contains the core functiions required to complete the simulations in the manuscript. The core model is found in _lhs_ode_gLV.m_ and is called within many of the other functions within the folder such as _change_parameter.m_, which runs the model for different parameter values.

Other functions in the folder help with assessment of analytical solutions of the model steady-states and use of the analytical solutions to explore the parameter space such as _sympbolic_solns.m_ and _calc_SS_stability.m_.


#### Global Sensitivity
The Global Sensitiivity folder has files created originally by the Kirschner Lab group as described in their 2008 paper, ["A Methodology For Performing Global Uncertainty And Sensitivity Analysis in Systems Biology"](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2570191/). Some minor modifications were made to help easily modify entries for the global sensitivity analysis that were specific to our model.

Examples of how to run the analysis are in the **Run_Files** folder.


