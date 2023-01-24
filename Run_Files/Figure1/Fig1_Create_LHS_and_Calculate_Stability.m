%% Fig1_run_global_SA.m
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Goal: Complete global sensitivity analysis given a set or possible
% parameter values using Latin Hypercube Sampling. The workspaces generated
% in this analysis are used for all downstream analyses.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Description of code workflow:
%   (1) Run the LHS-PRCC code
%           - Input: xlsx file
%           - Output: Model_LHS workspace
%   (2) Calculate SS stability using local stability analysis
%           - Input: Model_LHS workspace
%           - Output: SSConfig workspace
%                       + Figure of Model SS frequency
%   (3) Relate Model SS to CST Equilibrium Behavior (EB)
%           - Input SSConfig workspace
%           - Output: SSConfig workspace w/ CST EB appended
%                       + Figure of CST EB frequency
%
% REQUIRES: [Add folders to MATLAB path]
%   * Kirschner Lab Global Sensitivity Analysis Folder (Global_Sensitivity)
%   * ODE_scripts and General_scripts Folder 
%       * lhs_ode_gLV.m (ODE matlab function)
%       * Input excel file with LHS distributions
%       * analyze_global_SS.m
%       * analytical-base.mat
%       * calc_SS_stability.m
%           - requires "S" which is the models symbolic equations, "Jmat" which
%           is the Jacobian of the system of ODEs, "Type" which is the
%           formulation of ODEs (all of these values are in the
%           analytical-base.mat workspace)
%       * analyze_Global_CST_SS.m
%   * MATLAB TOOLBOXES: Parallel Computing Toolbox, Bioinformatics Toolbox,
%           Symbolic Math Toolbox
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee (chyylee@umich.edu)
% University of Michigan
% Jun 22, 2022
% Jan 23, 2023 (Updated to be more efficient)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% 1. Generate LHS Parameter Set Given an Excel File

flnm = 'lhs_settings_input.xlsx'; % Excel file with necessary parameters
[params,ICs] = pull_LHS_parameters(flnm);

NR = 500; % Number of "samples", our work uses 5000

paramMatrix = defineMatrix(params, NR,'parameter'); % Generates LHS parameter sets

%% 2. Get Predicted Model Steady-States & CST Classisfied Equilibrium Behavior

% LHS is stochastic, to use our generated parameter sets, please load the
% parameter space from "workspaces/SSConfig-Analysis-Model_LHS_10x.mat'".
% load('../workspaces/SSConfig-Analysis-Model_LHS_10x.mat','paramMatrix'),
% otherwise, you can run "analyze_Global_SS.m".
%   - Enter: paramMatrix (LHS generated parameter space: NR by number of
%   parameters double)
s
analyze_Global_SS(paramMatrix,{'NO','Li','oLB'}) % Call function that determines SS configurations











