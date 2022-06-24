%% Fig3_run_LHS_analysis_menses.m
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Goal: Assess the effect of perturbations that simulate menses. These
% simulations are based on the impact of biogenic amines associated with
% menses on growth rates and lactic acid production of different LB species
% in Borgorgna et al. (2021).
%
% The core code has the capibility to screen for many different
% combinations of parameter perturbations. Below, it is looking at how
% changes to 4 different parameters at varying degrees of parameter value
% alterations impacts community composition over time.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
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
%   * Running plot_temp_pert_results.m gives figure panels seen in Figure 3
%       with:
%           * Mean +/- 95% CI of resistant runs
%           * Mean +/- 95% CI of sensitive runs
%           * Volcano plot comparing resistant vs sensitive runs
%   
% REQUIRES:
%   * SSConfig-Analysis-Model_[].mat
%   * simulate_CST_EB_response.m
%       * change_parameter_men.m
%       * LHS_trace_visualize_new.m
%   * plot_temp_pert_results.m
%
% NOTE: Some GUI and command window prompts will appear. This is to help
% guide the user and ensure the correct inputs are selected.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% Jan 12, 2021
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% 1.) Loads required workspace and sets up base variables
fdr_loc = '../workspaces/'; 
load(strcat(fdr_loc,'SSConfig-Analysis-Model_LHS_10x.mat'), 'LHSmat','mat',...
    'StbleSS','S','Jmat','Type','SS_namesv','sp_names','param_names',...
    'all_nm')

% Prompt for CST Equilibrium type
[indx,~] = listdlg('PromptString',{'Pick Steady State Type:'},...
    'ListString',SS_namesv);
ss_id = indx;
[sel_idx] = find(all_nm == SS_namesv(indx));
sel_nets = LHSmat(sel_idx,:); %final parameter sets to use

%% 2. Run the perturbation analysis (Below is an example for Menses)
%  For menses simulations, select a LB species as the initial dominating
%  species. For BV therapy simulations, select NO/BV species as initial
%  dominating species.

% Prompts for dominate species (L. iners or oLB, depending on SS type)
[indx,tf] = listdlg('PromptString',{'Pick Initial Dominating Species:'},'ListString',sp_names);
sp_idx = indx; % idx of dominant species

% ###### CODE INPUTS ######
% Perturbation Length & Initial Condiitions
L = 7; % Duration of change
sp_p = 7; % Start time of modification
ep_p = sp_p + L; % End time of modification
time_post = 50; % Length of simulation after perturbation
dom_ab = 0.7; min_ab = (1 - dom_ab)/2; ybase = repmat(min_ab,[1,3]); % setting up initial conditions
ybase(sp_idx) = dom_ab; % setting up initial conditions

% Perturbation Methodology & Parameter Selection
perChange = "percentx"; % percentx (base + abs(base)*scalingFactor); plusx (base + scalingFactor); foldx (base*scalingFactor)
pidx = [2,3,7,10]; % index of growth rate parameters for L. iners and oLB and interactions terms Li/oLB -> BV
p1 = [-1 -0.5]; % ScalingFactor for pidx(1) - First index entered pidx [kgrow-Li]
p2 = p1; % ScalingFactor for pidx(2) [kgrow-oLB]
p3 = [0.5 1]; % ScalingFactor for pidx(3) [Li -> BV]
p4 = p3; % ScalingFactor for pidx4) [oLB -> BV]
vectorCell = {p1,p2,p3,p4}; % Save into cell (needed for code to run)

% Plotting & Naming
plotTraj = true; % Generate individual plots (can take a long time to run, enter false to speed up code)
stts = extractBetween(SS_namesv{ss_id},'[',']');
ss_type = string(strcat(extractBefore(SS_namesv{ss_id},':'),'-',join(stts,'-'))); % run name (SS configuration)
% ###### MODIFY RUN ABOVE ######

%% 3. Call the function
[dirName] = simulate_CST_EB_response(ybase,ss_type,StbleSS,sel_nets,S,Jmat,Type,...
    perChange,sp_p,ep_p,time_post,plotTraj,vectorCell,pidx);

%% 4. Plot the results for a desired perturbation
% GUI prompt will ask user for:
%   (1) Evaluation Time Point (often day 0 or day 30 after perturbation)
%   (2) Evaluation Threshold (60% - 0.6 relative abundance)
%   (3) Initial Dominating Species (LBinit or NOinit)
%   (4) Dose (will give a drop down menu for you to select)

plot_temp_pert_results(dirName)