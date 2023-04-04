%% G3_ABX_Run_1D_Bifurcation_and_Volcanos
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Goal: Assess how changing the growth rate of nAB in a manner similar to an
% ABX perturbation impacts the predicted equilibrium behavior. This can be
% extended to any parameter in the model.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% REQUIRED INPUTS:
%   * Model_LHS workspace (SSconfig-Analysis-HMP.mat')
%   * Folder and location information generated from
%       G3_Run_ABX_TimeSeries_Analysis.m
%
% REQUIRED FUNCTIONS
%   * SS_landscape_1D_v2.m 
%       * calc_SS_stability
%       * get_SS_info_3sp.m
%   * plot_1D_Bifurcation.m
%   * plot_Volcano.m
%
% OUTPUT:
%   * Workspace with 1D Bifurcation information for each EB type selected
%   * Figures for the 1D Bifurcation
%   * Figures with the volcano comparison on mono-stable and multi-stable
%       EB types
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% v1: Jan 22, 2022
% v2: Jan 21, 2023 (Updated to clean unneeded lines of code)
% v3: March 28, 2023 (converted code to no longer read in workspaces for
% each parameter change)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 1. Load Required Datan and Enter Bifurcation Parameters
fd_id = 1; % for manuscript analysis, there was only 1 folder generated

fdr_loc = '../G1_EB/';
load(strcat(fdr_loc,'SSConfig-Analysis-HMP.mat'),'S','Jmat','param_names')
fdr_names = {'ABXTimeSeries/nAB-HMP_-602lhs-day7-29-Mar-2023'};
ws_nm = '29-Mar-2023-1param-mod-0-for-7d-run.mat';
loc_name = strcat(fdr_names{fd_id},'/',ws_nm);

load(loc_name,'sel_nets','all_nm_CST')

N = 3; % number of speces
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

%% 2. Define Range for Bifurcation
% Select Bifurcation Parameter
pidx = 1; % growth of nAB

% Define Range and Number of Iterations (plusx)
pmin = 0; % minimal value
pmax = -5; % maximal value
pnum = 30; % number of iterations

%% 3. Pick Equilibrium behavior to analyze

SS_names_CST = unique(all_nm_CST);
[indx,tf] = listdlg('ListString',SS_names_CST);
sel_idx = [];
for i = 1:length(indx)
    tmp = all_nm_CST == SS_names_CST(indx(i)); % all_nm has been updated to all_nm_CST in some versions of code (Number of LHS samples x 1 string array)
    sel_idx = [sel_idx; find(tmp)];
end
fin_nets = sel_nets(sel_idx,:);
prepb = strrep(regexprep(SS_names_CST{indx},{':','[',']','/'},{'-','','-',''}),' ','');
fdr_nm = strcat(prepb,'-',date);


%% 4. Calculate the 1D Bifurcation

num_st = 1;
num_iter = size(fin_nets,1);
% Run the global code
SSdata = cell(num_iter,1);
parfor net_id = num_st:num_iter
    base_params = fin_nets(net_id,:);
    disp('~~~~~~~~~~~~~~~~~~')
    disp(strcat("RUNNING: ", " SET #", num2str(net_id)))
    disp('~~~~~~~~~~~~~~~~~~')
    [data_out,~,~] = SS_landscape_1D_v2(N,base_params,S,Jmat,...
        pidx,pnum,pmin,pmax);
    SSdata{net_id} = {data_out};

end

save(strcat(fdr_nm,'.mat'),'SSdata','pidx','pnum',...
    'pmin','pmax')
%% 5. Plot Bifurcation Result 
p_range = linspace(pmin,pmax,pnum);
plot_1D_Bifurcation_v2(SSdata,pnum,p_range,param_names,pidx)
hold on

xline(-2.64) % plot line where "abx" dose was located in perturbation analysis
xlim([-3,0]) % set limits
ylim([0,100]) % set limits

savefig(strcat(fdr_nm,'.fig'))

%% 6. Generate volcano plots

tit_list = {'volcano_1SS_nAB_vs_2SS_NO_or_Li';
    'volcano_1SS_nAB_vs_2SS_NO_or_oLB'};

comp_lists = [1,2; % comparison index: 1SS nAB to 2SS nAB or Li
    1,3];  % comparison index: 1SS nAB to 2SS nAB or oLB

SS_names_CST = unique(all_nm_CST);
for j = 1:length(tit_list)
    comp_list = comp_lists(j,:);
    tit = tit_list{j};
    sel_idx = {};
    for i = 1:length(comp_list)
        tmp = all_nm_CST == SS_names_CST(comp_lists(j,i)); % all_nm has been updated to all_nm_CST in some versions of code (Number of LHS samples x 1 string array)
        sel_idx{end+1} = find(tmp);
    end
    sel_nets1 = sel_nets(sel_idx{1},:);
    sel_nets2 = sel_nets(sel_idx{2},:);
    
    alpha = 0.01;
    offset = 0.1;
    [x,y]  = plot_Volcano(sel_nets1,sel_nets2,alpha,offset,param_names,SS_names_CST(comp_list));
   
    savefig(strcat(tit,'.fig'))
    close
end