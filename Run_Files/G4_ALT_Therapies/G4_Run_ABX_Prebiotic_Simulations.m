%% G4_Run_ABX_Prebiotic_Simulations.m
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Goal: Show how a combination of antibiotic (kgrow-NO) and prebiotic
% (affecting oLB, kgrow-oLB) impact treatment of nAB
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% INPUT:
%   * SSConfig workspace generated in /G1_EB/
%   * Data sets used in G3_ABX/ TimeSeries
%   * A CST Equilibrium Behavior of interest (will pull parameter sets of a
%       certain steady-state configuation type)
%   * Several parameters to run perturbation (see section 2)
%
% OUTPUT:
%   * Folder with:
%       * workspace with simulation results
%
% REQUIRES:
%   * SSConfig-Analysis-Model_[].mat
%   * simulate_CST_EB_response_displace.m
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% v1: Jun 22, 2022
% v2: Jan 22, 2023 (updated plotting function and call)
% v3: March 30, 2023 (streamlines post-processing of results)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 0. Loads required workspace and sets up base variables
fd_id = 1;

fdr_loc = '../G1_EB/';
load(strcat(fdr_loc,'SSConfig-Analysis-HMP.mat'),'S','Jmat','param_names')
fdr_names = {'../G3_ABX/ABXTimeSeries/nAB-HMP_-602lhs-day7-29-Mar-2023'};
ws_nm = '29-Mar-2023-1param-mod-0-for-7d-run.mat';
loc_name = strcat(fdr_names{fd_id},'/',ws_nm);
N = 3; % number of speces

load(loc_name,'sel_nets','all_nm_CST')
if ~exist('all_nm_CST','var')
    % Get EB Type Assignment
    all_nm_CST = [];
    for idx = 1:size(sel_nets,1)
        [run_mat,~,~,~] = calc_SS_stability(N,sel_nets(idx,:),S,Jmat);
        [nmf] = get_VALENCIA_class(run_mat); 
        all_nm_CST = [all_nm_CST; nmf];
    end
    save(loc_name,'all_nm_CST',"-append")
end

%% 1. Pull EB of Interest
SS_names_CST = unique(all_nm_CST);
[indx,~] = listdlg('ListString',SS_names_CST);
sel_idx = [];
for i = 1:length(indx)
    tmp = all_nm_CST == SS_names_CST(indx(i)); % all_nm has been updated to all_nm_CST in some versions of code (Number of LHS samples x 1 string array)
    sel_idx = [sel_idx; find(tmp)];
end
fin_nets = sel_nets(sel_idx,:);

% names output folder
prepb = strrep(regexprep(SS_names_CST{indx},{':','[',']','/'},{'-','','-',''}),' ','');
fdr_nm = strcat(prepb,'-',date);
ss_type = strcat('PreAbx-',fdr_nm); 
%% 2. Run the perturbation analysis (Below is an example for Menses)

sp_names = {'nAB','Li','oLB'};
% Prompts for dominate species (NO, depending on SS type)
[indx,tf] = listdlg('PromptString',{'Pick Initial Dominating Species:'},...
    'ListString',sp_names);
sp_idx = indx; % idx of dominant species

% ###### MODIFY BELOW ######
% Select parameters and values
perChange = "plusx"; % percentx (base + abs(base)*scalingFactor); plusx (base + scalingFactor); foldx (base*scalingFactor)
pidx = [1,3]; % index of growth rate parameters for L. iners and oLB and interactions terms Li/oLB -> BV
p1 = [0 -0.125 -0.25 -0.5 -1 -2 -2.64]; % Fold-change in growth rate Li
p2 = [0 0.125 0.25 0.5 1 2 2.64]; % Fold-change in growth rate Li

vectorCell = {p1,p2}; % Get combinations

% Determine simulation and perturbation length
L = 7; % Duration of change
sp_p = 7; % Start time of modification
ep_p = sp_p + L; % End time of modification
time_post = 60; % Length of simulation after perturbation

% Simulation initial conditions
dom_ab = -0.10; % Initial relative abundance of dominant species
ybase = [0.10 0.05 0.05];
ybase(sp_idx) = dom_ab;
plotTraj = false; % Generate individual plots (can take a long time to run)

% ###### MODIFY RUN ABOVE ######

%% 3. Call Simulation Function
[fldrnm] = simulate_CST_EB_response_displace(ybase,ss_type,fin_nets,S,Jmat,...
    perChange,sp_p,ep_p,time_post,plotTraj,vectorCell,pidx);


    