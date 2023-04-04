%% G2_MENSES_Run_TimeSeries_Analysis.m
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Goal: Simulate menses perturbations of varying strengths (no
% impact/control, light, moderate, and strong) across populations of women
% who are typically oLB or Li dominant.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% REQUIRES:
%   * SSConfig-Analysis_yourfilename.mat (-HMP.mat)
%   * remove_outlier_states.m
%   * simulate_CST_EB_response_displace.m
%       * change_parameter_men.m
%       * LHS_trace_visualize_new.m
%
% NOTE: Some GUI and command window prompts will appear. This is to help
% guide the user and ensure the correct inputs are selected.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% v1: Jan 12, 2021
% v2: Jan 21, 2023 (streamlined some code blocks)
% v3: Mar 29, 2023 (Pulls data matched to cliniclall EB Frequencies)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% 1.) Loads required workspace and sets up base variables
clear;clc;

ss_type = 'Li-States-HMP'; % run name (update based on what is run)
fdr_loc = '../G1_EB/'; 
load(strcat(fdr_loc,'SSConfig-Analysis-HMP.mat'))

% checks compatibility with previous naming conventions
if ~exist('SS_names_CST','var')
    SS_names_CST = SS_namesv;
    all_nm_CST = all_nm;
end

% Prompt for EB subtype and Extract Samples with EB
[indx,~] = listdlg('PromptString',{'Pick Steady State Type:'},...
    'ListString',SS_names_CST);
ss_id = indx;
sel_idx = [];
for i = 1:length(indx)
    [sel_idx] = [sel_idx;find(all_nm_CST == SS_names_CST(indx(i)))]; % all_nm has been updated to all_nm_CST in some versions of code (Number of LHS samples x 1 string array)
end
sel_nets = LHSmat(sel_idx,:); %final parameter sets to use

%% 2. Run the perturbation analysis 

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

vectorCell = {[0 0 0 0; % Control (No menses strength)
    -0.5 -0.5 0.5 0.5; % Light (Weak menses strength)
   -1 -1 1 1; % Moderate (Moderate menses strength)
   -2 -2 1 1]}; % Strong (Strong menses strength)

% Plotting & Naming
plotTraj = false; % Generate individual plots (can take a long time to run, so false skips step)

% ###### MODIFY RUN ABOVE ######
dom_ab = -0.10; % Percent displacement from SS
ybase = [0.10 0.05 0.05];
ybase(sp_idx) = dom_ab;


%% 3. Remove samples that undergo composition shifts without perturbation
tspan = [0 200];
rm_indx = remove_outlier_states(sel_nets,ybase,S,Jmat,tspan);

rmsum = sum(rm_indx);
disp(strcat(num2str(rmsum)," Sets Removed"))
%% 4. Run the perturbation
[dirName] = simulate_CST_EB_response_displace(ybase,ss_type,sel_nets(~rm_indx,:),S,Jmat,...
    perChange,sp_p,ep_p,time_post,plotTraj,vectorCell,pidx);



