%% run_LHS_analysis_template_menses.m
%
% Example code for the analysis in figures 3 and 4.
%
% Example of how to run the perturbation simulations, using menses as an
% example. The core code has the capibility to screen for many different
% combinations of parameter perturbations. Below, it is looking at how
% changes to 4 different parameters at varying degrees of parameter value
% alterations impacts community composition over time.
%
% REQUIRES:
%   * SSConfig-Analysis_yourfilename.mat
%   * simulate_CST_EB_response.m
%       * change_parameter_men.m
%       * LHS_trace_visualize_new.m
%   * check_comp (at end of scirpt)
%
% NOTE: Some GUI and command window prompts will appear. This is to help
% guide the user and ensure the correct inputs are selected.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% Jan 12, 2021
% Update: Jan 21, 2023 (streamlined some code blocks)
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

% Prompts for dominate species (L. iners or oLB, depending on EB type)
[indx,tf] = listdlg('PromptString',{'Pick Initial Dominating Species:'},'ListString',sp_names);
sp_idx = indx; % idx of dominant species

% ###### MODIFY BELOW ######
% Perturbation Length & Initial Condiitions
L = 7; % Duration of change
sp_p = 7; % Start time of modification
ep_p = sp_p + L; % End time of modification
time_post = 50; % Length of simulation after perturbation

% Perturbation Methodology & Parameter Selection
perChange = "plusfoldx"; %  plusfoldx (base + abs(base)*scalingFactor); plusx (base + scalingFactor); foldx (base*scalingFactor)
pidx = [2,3,7,10]; % index of growth rate parameters for L. iners and oLB and interactions terms Li/oLB -> NO
p1 = [-2 -0.5]; % perChange growth rate Li
p2 = p1; % perChange in growth rate oLB
p4 = [0.5 2]; % perChange in interaction of Li -> NO
p5 = p4; % perChange of interation of oLB -> NO
vectorCell = {p1,p2,p4,p5}; % Get combinations
dom_ab = 0.7; % Initial relative abundance of dominant species

% Plotting & Naming
plotTraj = false; % Generate individual plots (can take a long time to run, so false skips step)

% ###### MODIFY RUN ABOVE ######
min_ab = (1 - dom_ab)/2; ybase = repmat(min_ab,[1,3]); % setting up initial conditions
ybase(sp_idx) = dom_ab; % setting up initial conditions
stts = extractBetween(SS_names_CST{ss_id},'[',']');
ss_type = string(strcat(extractBefore(SS_names_CST{ss_id},':'),'-',join(stts,'-'))); % run name (SS configuration)

%% 3. Call the function
[dirName] = simulate_CST_EB_response(ybase,ss_type,sel_nets,S,Jmat,...
    perChange,sp_p,ep_p,time_post,plotTraj,vectorCell,pidx);

%% 4. plot

[Num_Sensitive,Num_Resilient,Vol_Prism,...
    Vol_Stats] = plot_Sensitive_vs_Resilient(dirName);