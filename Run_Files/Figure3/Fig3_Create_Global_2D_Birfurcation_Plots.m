%% Fig2_Create_Bifurcation_Plots.m
%
% GOAL: Of parameter sets with the same equilibrium behavior, analyze what
% the most common effect of a permenant perturbation has on the system
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% REQUIRED INPUTS:
%   * Model_LHS workspace
%   * Desired CST Equilibrium Behavior (GUI prompt)
%   * Desired model parameters to change
%   * Desired range of parameter changes (fold addition)
%   * Desired resolution (number of parameter combinations)
%
% REQUIRED FUNCTIONS
%   * SS_landscape_global_loop.m 
%       * calc_SS_stability
%   * plot_2D_Bifurcation.m
%
% OUTPUTS:
%   * Folder with a workspace saved for each parameter set
%   * Plots of most frequent equilbrium state, the frequency and the %BV
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% Jan 12, 2021
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% 1. GENERATE OR LOAD REQUIRED DATA

clear;clc;
fdr_loc = '../Figure1/';
load(strcat(fdr_loc,'SSConfig-Analysis-Model_LHS.mat'),'LHSmat',...
     'S','Jmat','SS_names_CST','param_names','all_nm','StbleSS')

N = 3; % number of speces

% Pick you parameters
p1 = [2,3]; % growth of iners, growth of oLB
p2 = [7,10];  % Li -> NO, oLB -> NO
    
% Value ranges (fold addition: baseline + p1*baseline)
p1min = -2.5;
p1max = 2.5;
p2min = 2.5;
p2max = -2.5;

% Number of parameter values to try (pnum values between p1min and p1max)
pnum = 11;

%% 2. Pull the Parameter Desired CST Equilibrium Behavvior

[indx,tf] = listdlg('ListString',SS_names_CST);
sel_idx = [];
for i = 1:length(indx)
    tmp = all_nm == SS_names_CST(indx(i));
    sel_idx = [sel_idx; find(tmp)];
end
sel_nets = LHSmat(sel_idx,:);
prepb = strrep(regexprep(SS_names_CST{indx},{':','[',']','/'},{'-','','-',''}),' ','');
fdr_nm = strcat(prepb,'-',date);
mkdir(fdr_nm)

%% 3. Run Global Bifurcation Analaysis
% WARNING: This code can take hours, even days to run. To speed up runs,
% change pnum to a smaller value or run the code in batches (change num_st
% and num_fun)

num_st = 1; % Start index
num_fin = size(sel_nets,1); % End index
for net_id = num_st:num_fin
    base_params = sel_nets(net_id,:);
    disp('~~~~~~~~~~~~~~~~~~')
    disp(strcat("RUNNING: ", " SET #", num2str(net_id)))
    disp('~~~~~~~~~~~~~~~~~~')
    [data_out,svnm] = SS_landscape_2D(N,base_params,param_names,S,Jmat,...
        p1,p2,pnum,p1min,p1max,p2min,p2max,net_id,fdr_nm);
end
disp('Done!')


%% 3. Plot results of bifurcation
% Enter desired folder name to plot
% fdr_nm ="1SS-Li-CST-III-20-Jan-2023/"; % Example folder name
plot_2D_Bifurcation(fdr_nm) % call plotting function

