%% Fig1_run_global_SA.m
%
% Generates workspace for downstream analysis.This can be completed in
% two parts:
%   (1) Run the LHS-PRCC code
%           - Input: xlsx file
%           - Output: workspace
%   (2) Calculate SS stability using local stability analysis
%
% REQUIRES:
%   * Kirschner Lab Global Sensitivity Analysis Folder
%   * lhs_ode_gLV.m (ODE matlab function)
%   * Input excel file with LHS distributions
%   * analyze_global_SS.m
%   * analytical-base.mat
%   * calc_SS_stability.m
%       - requires "S" which is the models symbolic equations, "Jmat" which
%       is the Jacobian of the system of ODEs, "Type" which is the
%       formulation of ODEs (all of these values are in the
%       analytical-base.mat workspace)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% Jan 12, 2021
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% 1. Run LHS-PRCC 
% This calls the Kirschner Group's LHS-PRCC code, but allows us to tell the
% code what settings we want to use and feed in the name of the xlsx file
% with our parmeter ranges.
% To use the LHS runs generate in our manuscript, skip to section #2

% @lhs_ode_gLV (ODE script/formulation used here)
lhs_ode_run_new('lhs_ode_settings_GUI')

%% 2. Get Predicted Steady-States
% 'Model_LHS_5000.mat' is the output generate from section #1 in our
%   (LHS code is stochastic, will not always generate exact same parameter
%   sets)

fdr_loc = '../workspaces'; % folder location
lhs_nm = 'Model_LHS_5000.mat'; % load worspace generated in #1
analyze_Global_SS(strcat(fdr_loc,'/',lhs_nm)) % Call function that determines SS configurations





