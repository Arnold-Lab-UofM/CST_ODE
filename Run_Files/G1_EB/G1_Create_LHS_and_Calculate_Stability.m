%% G1_Create_LHS_and_Calculate_Stability.m
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
%                       + Figure of CST EB frequency
%
% REQUIRES: [Add folders to MATLAB path]
%       * Input excel file with LHS distributions
%       * analyze_global_SS.m
%       * calc_SS_stability.m
%           - requires "S" which is the models symbolic equations, "Jmat" which
%           is the Jacobian of the system of ODEs
%       * analyze_Global_CST_SS.m
%   * MATLAB TOOLBOXES: Parallel Computing Toolbox, Bioinformatics Toolbox,
%           Symbolic Math Toolbox
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee (chyylee@umich.edu)
% University of Michigan
% Jun 22, 2022
% Jan 23, 2023 (Updated to be more efficient)
% Apr 4, 2023 (Updated to be more reproducible)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% 1. Generate LHS Parameter Set Given an Excel File
rng(1);
flnm = 'lhs_settings_input_Litest.xlsx'; % Excel file with necessary parameters
[params,ICs] = pull_LHS_parameters(flnm);

NR = 5000; % Number of "samples", our work uses 5000

paramMatrix = defineMatrix(params, NR,'parameter'); % Generates LHS parameter sets

%% 2. Get Predicted Model Steady-States & CST Classisfied Equilibrium Behavior

% LHS is stochastic, to use our generated parameter sets, please load the
% parameter space from "workspaces/SSConfig-Analysis-Model_LHS_10x.mat'".
% load('../workspaces/SSConfig-Analysis-Model_LHS_10x.mat','paramMatrix'),
% otherwise, you can run "analyze_Global_SS.m".
%   - Enter: paramMatrix (LHS generated parameter space: NR by number of
%   parameters double)
%
% Note: Older versions of code may have 'NO' rather than 'nAB'

% load('../workspaces/SSConfig-Analysis-Model_LHS_10x.mat','paramMatrix')
analyze_Global_SS(paramMatrix,{'nAB','Li','oLB'}) % Call function that determines SS configurations






