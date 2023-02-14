%% Fig5_run_LHS_analysis_BV_abx_dose_dur.m
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Goal: Simulate an antibiotic therapy (metronidazole) dose and duration
% by impacting the growth rate of the NO population and length of
% perturbation.
%
% Example code for the analysis in figure 5B.
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
%   * plot_dosedur.m takes results and gives heatmap of results 
%
% REQUIRES:
%   * SSConfig-Analysis-Model_[].mat
%   * simulate_CST_EB_response.m
%       * change_parameter_men.m
%       * LHS_trace_visualize_new.m
%   * plot_dosedur.m
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% Jan 12, 2021
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 1. Loads required workspace and sets up base variables
clear;clc;
fdr_loc = '../Figure1/'; 
load(strcat(fdr_loc,'SSConfig-Analysis-Model_LHS.mat'), 'LHSmat','StbleSS','S',...
    'Jmat','SS_names_CST','sp_names','param_names','all_nm_CST')

% Prompt for CST Equilibrium Behavior:
[indx,tf] = listdlg('PromptString',{'Pick Steady State Type:'},...
    'ListString',SS_names_CST);
ss_id = indx;
[sel_idx] = find(all_nm_CST == SS_names_CST(indx)); % all_nm has been updated to all_nm_CST in some versions of code (Number of LHS samples x 1 string array)
sel_nets = LHSmat(sel_idx,:); %final parameter sets to use
%% 2. Run Dose and Duration Simulations 

% Input Code Input Parameters
pidx = [1]; % index of parameter of interest
Lin = [1 7 14 21 60]; % duration (days)
var_vals = [-3 -2 -1 -0.75 -0.5 0]; % dose
time_post = 50; % how far after alteration to analyze
perChange = "plusx"; % (same as percentx: base + abs(base)*scalingFactor)

% Initial Conditions
indx = 1; sp_idx = indx; % idx of dominant species
dom_ab = 0.7; min_ab = (1 - dom_ab)/2; ybase = repmat(min_ab,[1,3]);
ybase(sp_idx) = dom_ab;

% Other Conditions
plotTraj = false; % Generate individual plots (can take a long time to run)
ss_type = strrep(SS_names_CST{ss_id},'/','');

% Loop to look at combinations
for v = 1:length(var_vals)
    p1 = [var_vals(v)]; % Change in growth rate 1
    vectorCell = {p1}; % Get combinations
    for i = 1:length(Lin)
        L = Lin(i);
        sp_p = 7; % Start time of modification
        ep_p = sp_p + L; % End time of modification
%         [arm] = simulate_LHS_time_response(ybase,ss_type,...
%             sel_nets,S,Jmat,perChange,sp_p,ep_p,time_post,vectorCell,pidx);
        [arm] = simulate_CST_EB_response(ybase,ss_type,sel_nets,S,Jmat,perChange,sp_p,ep_p,time_post,plotTraj,vectorCell,pidx);
    end
end
%% 3. Plot Results
dose = var_vals; % selected doses
duration = Lin; % selected duration
ev_point = 0; % Time point for evaluation (post regimen)
fldrnm = '1SS: [NO] CST-IV_1D_k_grow-NO_23-Jan-2023'; % Example folder

% run plot file
[percent_switch,total_runs] = plot_DoseRegimens_Heatmap(fldrnm,dose,duration,ev_point);
