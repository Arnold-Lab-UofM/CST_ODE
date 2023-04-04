%% G4_Run_Dose_Duration_Simulations.m
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Goal: Simulate an antibiotic therapy (metronidazole) dose and duration
% by impacting the growth rate of the NO population and length of
% perturbation.
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
%       * change_parameter.m
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% v1: Jan 12, 2021
% v2: Apr 4, 2023 (Updated with new initial condition methodology)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 0. Loads required workspace and sets up base variables
fd_id = 1;

fdr_loc = '../G1_EB/';
load(strcat(fdr_loc,'SSConfig-Analysis-HMP.mat'),'S','Jmat','param_names')
fdr_names = {'../G3_ABX/ABXTimeSeries/nAB-HMP_-602lhs-day7-29-Mar-2023'};
ws_nm = '29-Mar-2023-1param-mod-0-for-7d-run.mat';
loc_name = strcat(fdr_names{fd_id},'/',ws_nm);

load(loc_name,'sel_nets','all_nm_CST')
if ~exist('all_nm_CST','var')
    % Get EB Type Assignment
    N = 3; % number of species
    all_nm_CST = [];
    for idx = 1:size(sel_nets,1)
        [run_mat,~,~,~] = calc_SS_stability(N,sel_nets(idx,:),S,Jmat);
        [nmf] = get_VALENCIA_class(run_mat); 
        all_nm_CST = [all_nm_CST; nmf];
    end
    save(loc_name,'all_nm_CST',"-append")
end
%% 1. Indicate EB to Analyze
SS_names_CST = unique(all_nm_CST);
[indx,tf] = listdlg('ListString',SS_names_CST);
sel_idx = [];
for i = 1:length(indx)
    tmp = all_nm_CST == SS_names_CST(indx(i)); % all_nm has been updated to all_nm_CST in some versions of code (Number of LHS samples x 1 string array)
    sel_idx = [sel_idx; find(tmp)];
end
fin_nets = sel_nets(sel_idx,:);

% names output folder
prepb = strrep(regexprep(SS_names_CST{indx},{':','[',']','/'},{'-','','-',''}),' ','');
fdr_nm = strcat(prepb,'-',date);
ss_type = strcat('DoseDur-',fdr_nm); 

%% 2. Set Up Simulation Inputs
sp_idx = 1; % nAB start for BV

% Input Code Input Parameters
pidx = [1]; % index of parameter of interest
duration = [1 7 14 21 60]; % duration (days)
doses = [-2.64 -2 -1 -0.5 -0.25 -0.125 0]; % dose
time_post = 50; % how far after alteration to analyze
perChange = "plusx"; % (same as percentx: base + abs(base)*scalingFactor)

% Simulation initial conditions
dom_ab = -0.10; % Initial relative abundance of dominant species
ybase = [0.10 0.05 0.05];
ybase(sp_idx) = dom_ab;

% Other Conditions
plotTraj = false; % Generate individual plots (can take a long time to run)

%% 3. Run Simulation
vectorCell = {[doses']};
for i = 1:length(duration)
    L = duration(i);
    sp_p = 7; % Start time of modification
    ep_p = sp_p + L; % End time of modification
    [arm] = simulate_CST_EB_response_displace(ybase,ss_type,fin_nets,S,Jmat,perChange,sp_p,ep_p,time_post,plotTraj,vectorCell,pidx);
end


