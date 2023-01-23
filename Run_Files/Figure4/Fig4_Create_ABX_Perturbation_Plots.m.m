%% Fig4_run_abx_perturbation.m
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Goal: Simulate an antibiotic therapy (metronidazole) by impacting the
% growth rate of the NO population. The parameter changes were based on
% values calculated from Mayer et al. 2015, which analyzed changes in
% abudance of BV associated bacteria after an antibiotic treatment for BV.
%
% Example code for the analysis in figures 4.
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
%   * Running plot_temp_pert_results.m gives figure panels seen in Figure 3
%       with:
%           * Mean +/- 95% CI of resistant runs
%           * Mean +/- 95% CI of sensitive runs
%           * Volcano plot comparing resistant vs sensitive runs
%   
% REQUIRES:
%   * SSConfig-Analysis-Model_[].mat
%   * simulate_CST_EB_response.m
%       * change_parameter.m
%       * plot_CST_EB_response.m
%   * plot_Sensitive_vs_Resilient.m
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% v1: Jan 12, 2021
% v2: Jan 20, 2023 (Update to new naming conventions)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% 1.) Loads required workspace and sets up base variables
fdr_loc = '../Figure1/'; 
load(strcat(fdr_loc,'SSConfig-Analysis-Model_LHS.mat'), 'LHSmat','StbleSS','S',...
    'Jmat','SS_names_CST','sp_names','param_names','all_nm')

% Prompt for SS type:
[indx,~] = listdlg('PromptString',{'Pick Steady State Type:'},...
    'ListString',SS_names_CST);
ss_id = indx;
[sel_idx] = find(all_nm == SS_names_CST(indx));
sel_nets = LHSmat(sel_idx,:); %final parameter sets to use

%% 2. Run the perturbation analysis (Below is an example for Menses)
%   To modify for ABX, select BV as the initial dominating species, change
%   perChange to plusx (now will subtract rather than fold change), change
%   pidx = 1 (kgrow_BV), The average value across 12 patients and 6
%   genuses of BVAB is 2.64 h-1 in Mayer et al. 2015

% Prompts for dominate species (NO/BV for ABX simulations)
[indx,tf] = listdlg('PromptString',{'Pick Initial Dominating Species:'},...
    'ListString',sp_names);
sp_idx = indx; % idx of dominant species

% ###### CODE INPUTS ######
% Perturbation Length & Initial Condiitions
L = 7; % Duration of change
sp_p = 7; % Start time of modification
ep_p = sp_p + L; % End time of modification
time_post = 60; % Length of simulation after perturbation
dom_ab = 0.7; min_ab = (1 - dom_ab)/2; % setting up initial conditions
ybase = repmat(min_ab,[1,3]); ybase(sp_idx) = dom_ab; % setting up initial conditions

% Perturbation Methodology & Parameter Selection
perChange = "plusx"; % percentx (base + abs(base)*scalingFactor); plusx (base + scalingFactor); foldx (base*scalingFactor)
pidx = [1]; % index of growth rate parameters for NO
p1 = [0  -0.5 -1 -2 -2.64 -5]; % Scaling Factor for NO growth (death)
vectorCell = {p1}; % Input as cell

% Plotting & Naming
plotTraj = false; % Generate individual plots (can take a long time to run)
stts = extractBetween(SS_names_CST{ss_id},'[',']');
ss_type = string(strcat(extractBefore(SS_names_CST{ss_id},':'),'-',join(stts,'-'))); % run name (SS configuration)
% ###### MODIFY RUN ABOVE ######

%% 3. Call the function to run temp. perturbations
[fdr_nm] = simulate_CST_EB_response(ybase,ss_type,sel_nets,S,Jmat,...
    perChange,sp_p,ep_p,time_post,plotTraj,vectorCell,pidx);

%% 4. Plot results
% fdr_nm = '2SS-NO-Li_1D_k_grow-BV_22-Jun-2022'; % enter folder name
% GUI prompt will ask user for:
%   (1) Evaluation Time Point (often day 0 or day 30 after perturbation)
%   (2) Evaluation Threshold (60% - 0.6 relative abundance)
%   (3) Initial Dominating Species (LBinit or NOinit)
%   (4) Dose (will give a drop down menu for you to select)

[Num_Sensitive,Num_Resilient,Vol_Prism,...
    Vol_Stats] = plot_Sensitive_vs_Resilient(fdr_nm);
