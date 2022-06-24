%% Fig1_run_global_SA.m
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Goal: Complete global sensitivity analysis given a set or possible
% parameter values using Latin Hypercube Sampling. The workspaces generated
% in this analysis are used for all downstream analyses.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Description of code workflow:
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
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee (chyylee@umich.edu)
% University of Michigan
% Jun 22, 2022
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% 1. Run LHS-PRCC 
% This calls the Kirschner Group's LHS-PRCC code, but allows us to tell the
% code what settings we want to use and feed in the name of the xlsx file
% with our parmeter ranges.
% To use the LHS runs generate in our manuscript, skip to section #2

% @lhs_ode_gLV (ODE script/formulation used here)
lhs_ode_run_new('lhs_ode_settings_GUI')

%% 2. Get Predicted Model Steady-States
% 'Model_LHS.mat' is the output generate from section #1 in our
%   (LHS code is stochastic, will not always generate exact same parameter
%   sets)

fdr_loc = ''; % folder location
lhs_nm = 'Model_LHS.mat'; % load Model_LHS workspace in #1 
analyze_Global_SS(strcat(fdr_loc,lhs_nm),false) % Call function that determines SS configurations

%% 3. Get CST Equilibrium Behavior
% To relate model SS to clinical behavior, a nearest centroid classifier
% bins each model SS to a CST from centroid published in France et al.,
% 2020 (VALENCIA)

fn = 'SSConfig-Analysis-Model_LHS.mat';
analyze_Global_CST_SS(fn)




