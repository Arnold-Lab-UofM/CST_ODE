%% Fig4_Create_1D_Global_Bifurcation_Plots.m
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% REQUIRED INPUTS:
%   * Model_LHS workspace
%
% REQUIRED FUNCTIONS
%   * SS_landscape_1D.m 
%       * calc_SS_stability
%       * get_SS_info_3sp.m
%   * Plot_1D_Global_Bifurcation.m
%
% OUTPUT:
%   * Folder with result of each parameter change
%   * Figure with results (frequency of each CST equilibrium behavior)
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% Jan 12, 2021
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%% 1. GENERATE OR LOAD REQUIRED DATA

clear;clc;
fdr_loc = '../workspaces';
load(strcat(fdr_loc,'/SSConfig-Analysis-Model_LHS_10x.mat'),'LHSmat',...
    'mat', 'S','Jmat','Type','colors', ...
    'param_names','all_nm','StbleSS')

% Parameter
p2 = [1]; % growth of BV

% Range and iterations
p2min = 0; 
p2max = -5;
pnum2 = 30;

%%
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

%% 2. CALL 2D BIFURCATION FUNCTION
num_st = 1;
num_iter = 20;

p1 = p2; 
p1min = 0;
p1max = 0.1;
pnum1 = 2;
for i = num_st:num_iter
    net_id = i;
    base_params = sel_nets(net_id,:);
    disp('~~~~~~~~~~~~~~~~~~')
    disp(strcat("RUNNING: ", " SET #", num2str(net_id)))
    disp('~~~~~~~~~~~~~~~~~~')

    [SS_map,data_out,sum_table,svnm] = SS_landscape_1D(3,base_params,param_names,S,Jmat,Type,...
        colors,indx,p1,p2,pnum1,pnum2,p1min,p1max,p2min,p2max,net_id,fdr_nm);
end

%% 3. Plot Result
Plot_1D_Global_Bifurcation(fdr_nm)