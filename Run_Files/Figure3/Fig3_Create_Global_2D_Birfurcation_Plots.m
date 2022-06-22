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
%       * get_SS_info_3sp.m
%   * Plot_Global_Bifurcation.m
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
fdr_loc = '../workspaces/';
load(strcat(fdr_loc,'SSConfig-Analysis-Model_LHS_10x.mat'),'LHSmat',...
    'mat', 'S','Jmat','Type','colors', ...
    'param_names','all_nm','StbleSS')

% Pick you paraameters
p1 = [2,3]; % growth of iners, growth of oLB
p2 = [7,10];  % Li -> NO, oLB -> NO
    
% Value ranges (fold addition: baseline + p1*baseline)
p1min = -2.5;
p1max = 2.5;
p2min = 2.5;
p2max = -2.5;

% Number of parameter values to try (51 values between p1min and p1maax)
pnum = 51;

%% 2. Pull the Parameter Desired CST Equilibrium Behavvior

CST_state = {'1SS: [Li] CST-III';'1SS: [oLB] CST-I/II/V';'1SS: [NO] CST-IV';'2SS: [NO] CST-IV or [oLB] CST-I/II/V';'2SS: [Li] CST-III or [oLB] CST-I/II/V';'2SS: [Li] CST-III or [Li] CST-III';'2SS: [NO] CST-IV or [Li] CST-III';'3SS: [NO] CST-IV or [Li] CST-III or [oLB] CST-I/II/V';'2SS: [oLB] CST-I/II/V or [oLB] CST-I/II/V';'2SS: [NO] CST-IV or [NO] CST-IV'};
[indx,tf] = listdlg('ListString',CST_state);
sel_idx = [];
for i = 1:length(indx)
    tmp = all_nm == CST_state(indx(i));
    sel_idx = [sel_idx; find(tmp)];
end
sel_nets = LHSmat(sel_idx,:);
CST_state(indx);
prepb = strrep(regexprep(CST_state{indx},{':','[',']','/'},{'-','','-',''}),' ','');
fdr_nm = strcat(prepb,'-',date);
mkdir(fdr_nm)

%% 3. Run Global Bifurcation Analaysis
% WARNING: This code can take hours, even days to run. To speed up runs,
% change pnum to a smaller value or run the code in batches (change num_st
% and num_fun)

num_st = 1; % Start index
num_fin = size(sel_nets,1); % End index


for i = num_st:num_fin
    net_id = i;
    base_params = sel_nets(net_id,:);
    disp('~~~~~~~~~~~~~~~~~~')
    disp(strcat("RUNNING: ", " SET #", num2str(net_id)))
    disp('~~~~~~~~~~~~~~~~~~')
    [SS_map,data_out,sum_table,svnm] = SS_landscape_global_loop(3,base_params,param_names,S,Jmat,Type,...
        colors,indx,p1,p2,pnum,p1min,p1max,p2min,p2max,net_id,fdr_nm);
end
disp('Done!')


%% 3. Plot results of bifurcation
% Enter desired folder name to plot

Plot_Global_Bifurcation(fdr_nm)

