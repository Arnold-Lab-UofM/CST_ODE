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
%       * change_parameter_men.m
%       * LHS_trace_visualize_new.m
%   * [Results,~,~] = plot_comboresults(fldrnm,ev_val)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% Jun 22, 2022
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 1.) Loads required workspace and sets up base variables
fdr_loc = '../workspaces/'; 
load(strcat(fdr_loc,'SSConfig-Analysis-Model_LHS_10x.mat'), 'LHSmat','mat','StbleSS','S',...
    'Jmat','Type','SS_namesv','sp_names','param_names','all_nm')

% Prompt for SS type:
[indx,~] = listdlg('PromptString',{'Pick Steady State Type:'},...
    'ListString',SS_namesv);
ss_id = indx;
[sel_idx] = find(all_nm == SS_namesv(indx));
sel_nets = LHSmat(sel_idx,:); %final parameter sets to use

%% 2. Run the perturbation analysis (Below is an example for Menses)
%   To modify for ABX, select BV as the initial dominating species, change
%   perChange to false (now will subtract rather than fold change), change
%   pidx = 1 (kgrow_BV), p1 = 0.95 d^-1 (change in growth rate), and
%   vectorCell should now only have one entry {p1}, make sure to change
%   "init_state" above to 'BV' The aveage value across 12 patients and 6
%   genuses of BVAB is 2.64 h-1

% Prompts for dominate species (L. iners or oLB, depending on SS type)
[indx,tf] = listdlg('PromptString',{'Pick Initial Dominating Species:'},'ListString',sp_names);
sp_idx = indx; % idx of dominant species

% ###### MODIFY BELOW ######
pidx = [1,3]; % index of growth rate parameters for L. iners and oLB and interactions terms Li/oLB -> BV
p1 = [0 -0.1 -0.25 -0.5 -1 -2 -2.5 -3 -5]; % Fold-change in growth rate Li
p2 = [0 0.1 0.25 0.5 1 2 2.5 3 5]; % Fold-change in growth rate Li
vectorCell = {p1,p2}; % Get combinations

L = 7; % Duration of change
sp_p = 7; % Start time of modification
ep_p = sp_p + L; % End time of modification
plotTraj = false; % Generate individual plots (can take a long time to run)
perChange = "percentx"; % percentx (base + abs(base)*scalingFactor); plusx (base + scalingFactor); foldx (base*scalingFactor)
time_post = 60; % Length of simulation after perturbation
dom_ab = 0.7; % setting up initial conditions
min_ab = (1 - dom_ab)/2; % setting up initial conditions
ybase = repmat(min_ab,[1,3]); % setting up initial conditions
ybase(sp_idx) = dom_ab; % setting up initial conditions
stts = extractBetween(SS_namesv{ss_id},'[',']');
ss_type = string(strcat(extractBefore(SS_namesv{ss_id},':'),'-',join(stts,'-'))); % run name (SS configuration)
% ###### MODIFY RUN ABOVE ######

% Call the function
[fldrnm] = simulate_CST_EB_response(ybase,ss_type,StbleSS,sel_nets,S,Jmat,Type,...
    perChange,sp_p,ep_p,time_post,plotTraj,vectorCell,pidx);

%% 3. Plot Results
% fldrnm = '2SS-NO-Li_2D_k_grow-BVk_grow-oLB_22-Jun-2022';
ev_val = 0; % days post regimen
[Results,~,~] = plot_comboresults(fldrnm,ev_val);
    