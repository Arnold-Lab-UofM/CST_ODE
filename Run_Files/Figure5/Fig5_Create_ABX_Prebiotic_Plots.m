%% run_LHS_analysis_abx_prebiotic_combo.m
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Goal: Show how a combination of antibiotic (kgrow-NO) and prebiotic
% (affecting oLB, kgrow-oLB) impact treatment of NO. This is similar to
% code in Figure 4, but now two variables are being varied at one time andd
% a heatmap of the treatment efficacy is returned for evaluated at a user
% specified time point.
%
% Example code for the analysis in figure 5A.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% INPUT:
%   * SSConfig workspace generated in Figure 1
%   * A CST Equilibrium Behavior of interest (will pull parameter sets of a
%       certain steady-state configuation type)
%   * Several parameters to run perturbation (see section 2)
%
% OUTPUT:
%   * Folder with:
%       * workspace with simulation results
%       * summary figure of results (evaluated at run end point)
%   * Running plot_comboresults.m gives heatmap of results
%
% REQUIRES:
%   * SSConfig-Analysis-Model_[].mat
%   * simulate_CST_EB_response.m
%       * change_parameter.m
%       * plot_CST_EB_response.m
%   * [Results,~,~] = plot_Combination_Heatmap(fldrnm,ev_val)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% v1: Jun 22, 2022
% v2: Jan 22, 2023 (updated plotting function and call)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 1.) Loads required workspace and sets up base variables
fdr_loc = '../Figure1/'; 
load(strcat(fdr_loc,'SSConfig-Analysis-Model_LHS.mat'), 'LHSmat','StbleSS','S',...
    'Jmat','SS_names_CST','sp_names','param_names','all_nm_CST')

% Prompt for SS type:
[indx,~] = listdlg('PromptString',{'Pick Steady State Type:'},...
    'ListString',SS_names_CST);
ss_id = indx;
[sel_idx] = find(all_nm_CST == SS_names_CST(indx)); % all_nm has been updated to all_nm_CST in some versions of code (Number of LHS samples x 1 string array)
sel_nets = LHSmat(sel_idx,:); %final parameter sets to use

%% 2. Run the perturbation analysis (Below is an example for Menses)

% Prompts for dominate species (NO, depending on SS type)
[indx,tf] = listdlg('PromptString',{'Pick Initial Dominating Species:'},...
    'ListString',sp_names);
sp_idx = indx; % idx of dominant species

% ###### MODIFY BELOW ######
% Select parameters and values
perChange = "plusfoldx"; % percentx (base + abs(base)*scalingFactor); plusx (base + scalingFactor); foldx (base*scalingFactor)
pidx = [1,3]; % index of growth rate parameters for L. iners and oLB and interactions terms Li/oLB -> BV
p1 = [0 -0.1 -0.25 -0.5 -1 -2 -2.5 -3 -5]; % Fold-change in growth rate Li
p2 = [0 0.1 0.25 0.5 1 2 2.5 3 5]; % Fold-change in growth rate Li
vectorCell = {p1,p2}; % Get combinations

% Determine simulation and perturbation length
L = 7; % Duration of change
sp_p = 7; % Start time of modification
ep_p = sp_p + L; % End time of modification
time_post = 60; % Length of simulation after perturbation

% Simulation initial conditions
dom_ab = 0.7; % setting up initial conditions
min_ab = (1 - dom_ab)/2; % setting up initial conditions
ybase = repmat(min_ab,[1,3]); % setting up initial conditions
ybase(sp_idx) = dom_ab; % setting up initial conditions

plotTraj = false; % Generate individual plots (can take a long time to run)
stts = extractBetween(SS_names_CST{ss_id},'[',']');
ss_type = string(strcat(extractBefore(SS_names_CST{ss_id},':'),'-',...
    join(stts,'-'))); % run name (SS configuration)
% ###### MODIFY RUN ABOVE ######

% Call the function
[fldrnm] = simulate_CST_EB_response(ybase,ss_type,sel_nets,S,Jmat,...
    perChange,sp_p,ep_p,time_post,plotTraj,vectorCell,pidx);

%% 3. Plot Results
ev_val = 0; % days post regimen
[results,total_runs,param_values1,...
    param_values2] = plot_Combination_Heatmap(fldrnm,ev_val);
    