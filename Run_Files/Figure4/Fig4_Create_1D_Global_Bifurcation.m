%% Fig4_run_1D_Global_Bifurcation_Plots.m
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Goal: Assess how changing the growth rate of NO in a manner similar to an
% ABX perturbation impacts the predicted equilibrium behavior. This can be
% extended to any parameter in the model.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% REQUIRED INPUTS:
%   * Model_LHS workspace
%
% REQUIRED FUNCTIONS
%   * SS_landscape_1D.m 
%       * calc_SS_stability
%       * get_SS_info_3sp.m
%   * plot_1D_Bifurcation.m
%
% OUTPUT:
%   * Folder with result of each parameter change
%   * Figure with results (frequency of each CST equilibrium behavior)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% v1: Jan 22, 2022
% v2: Jan 21, 2023 (Updated to clean unneeded lines of code)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 1. Load Required Datan and Enter Bifurcation Parameters
fdr_loc = '../Figure1';
load(strcat(fdr_loc,'/SSConfig-Analysis-Model_LHS.mat'),'LHSmat',...
    'S','Jmat', 'param_names','all_nm_CST','StbleSS','SS_names_CST')

% ###### Modify here #######
% Select Bifurcation Parameter
pidx = 1; % growth of BV

% Define Range and Number of Iterations (plusx)
pmin = 0; % minimal value
pmax = -5; % maximal value
pnum = 30; % number of iterations
% ###### Modify here #######
%% 2. Pick Equilibrium behavior to analyze

[indx,tf] = listdlg('ListString',SS_names_CST);
sel_idx = [];
for i = 1:length(indx)
    tmp = all_nm_CST == SS_names_CST(indx(i)); % all_nm has been updated to all_nm_CST in some versions of code (Number of LHS samples x 1 string array)
    sel_idx = [sel_idx; find(tmp)];
end
sel_nets = LHSmat(sel_idx,:);
prepb = strrep(regexprep(SS_names_CST{indx},{':','[',']','/'},{'-','','-',''}),' ','');
fdr_nm = strcat(prepb,'-',date);
mkdir(fdr_nm)

%% 2. Call 1D Bifurcation Code

num_st = 1;
num_iter = size(sel_nets,1);
N = 3; % number of species
% Run the global code
for net_id = num_st:num_iter
    base_params = sel_nets(net_id,:);
    disp('~~~~~~~~~~~~~~~~~~')
    disp(strcat("RUNNING: ", " SET #", num2str(net_id)))
    disp('~~~~~~~~~~~~~~~~~~')
    [data_out,svnm] = SS_landscape_1D(N,base_params,param_names,S,Jmat,...
        pidx,pnum,pmin,pmax,net_id,fdr_nm);
end
%% 3. Plot Result (require folder name of results to run "fdr_nm")
plot_1D_Bifurcation(fdr_nm)
hold on

% Add Figure 4 specific line: "abx" dose
xline(-2.64) % plot line where "abx" dose was located
xlim([-3,0]) % set limits
ylim([0,100]) % set limits