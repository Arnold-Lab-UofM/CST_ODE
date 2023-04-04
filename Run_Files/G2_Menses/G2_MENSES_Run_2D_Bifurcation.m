%% G2_MENSES_Run_2D_Bifurcation.m
%
% GOAL: Of parameter sets with the same equilibrium behavior, analyze what
% the most common effect of a permenant perturbation has on the system
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% REQUIRED INPUTS:
%   * Model_LHS workspaces
%       - Uses workspaces generated from
%       G2_MENSES_Run_TimeSeries_Analyses.m (will need to indicate names
%       and locations of folders)
%   * Desired CST Equilibrium Behavior (GUI prompt)
%   * Desired model parameters to change
%   * Desired range of parameter changes (fold addition)
%   * Desired resolution (number of parameter combinations, pnum)
%
% REQUIRED FUNCTIONS
%   * plot_2D_Bifurcation.m
%
% OUTPUTS:
%   * Folder with a workspace saved for each parameter set
%   * Plots of most frequent equilbrium state, the frequency and the %BV
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% v1: Jan 12, 2021
% v2: Apr 4, 2023 (uses in silico HMP population)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% 0. GENERATE OR LOAD REQUIRED DATA

clear;clc;
fd_id = 1; % 1 = oLB states, 2 = Li states

fdr_loc = '../G1_EB/';
load(strcat(fdr_loc,'SSConfig-Analysis-HMP.mat'),'S','Jmat','param_names')

output_folder = 'MensesTimeSeries/';
fdr_names = {'oLB-States-HMP_-446lhs-day7-29-Mar-2023';
    'Li-States-HMP_-309lhs-day7-29-Mar-2023'};

ws_nm = '29-Mar-2023-4param-mod-0  0  0  0-for-7d-run.mat';
loc_name = strcat(output_folder,fdr_names{fd_id},'/',ws_nm);
load(loc_name,'sel_nets','all_nm_CST')

 % Get EB Type Assignment
 N = 3; % number of speces
if ~exist('all_nm_CST','var')
    all_nm_CST = [];
    for idx = 1:size(sel_nets,1)
        [run_mat,~,~,~] = calc_SS_stability(N,sel_nets(idx,:),S,Jmat);
        [nmf] = get_VALENCIA_class(run_mat); 
        all_nm_CST = [all_nm_CST; nmf];
    end
    save(loc_name,'all_nm_CST',"-append")
end
%% 1. INDICATE PARAMETER ALTERATIONS

% Pick you parameters
p1 = [2,3]; % growth of iners, growth of oLB
p2 = [7,10];  % Li -> NO, oLB -> NO
    
% Value ranges (fold addition: baseline + p1*baseline)
p1min = -2.5;
p1max = 2.5;
p2min = 2.5;
p2max = -2.5;

% Number of parameter values to try (pnum values between p1min and p1max)
pnum = 51;

%% 2. Pull the Parameter Desired CST Equilibrium Behavvior
SS_names_CST = unique(all_nm_CST);
[indx,tf] = listdlg('ListString',SS_names_CST);
sel_idx = [];
for i = 1:length(indx)
    tmp = all_nm_CST == SS_names_CST(indx(i)); % all_nm has been updated to all_nm_CST in some versions of code (Number of LHS samples x 1 string array)
    sel_idx = [sel_idx; find(tmp)];
end
fin_nets = sel_nets(sel_idx,:);

% output name
prepb = strrep(regexprep(SS_names_CST{indx},{':','[',']','/'},{'-','','-',''}),' ','');
fdr_nm = strcat(prepb,'-',date);
matname = strcat('x2D_Bif_',fdr_nm);

%% 3. Run Global Bifurcation Analaysis
% WARNING: This code can take hours, even days to run. To speed up runs,
% change pnum to a smaller value or run the code in batches (change num_st
% and num_fun)

num_st = 1; % Start index
num_fin = size(fin_nets,1); % End index
alldata_out = cell(num_fin,1);
c = 1;
for net_id = num_st:num_fin
    base_params = fin_nets(net_id,:);
    disp('~~~~~~~~~~~~~~~~~~')
    disp(strcat("RUNNING: ", " SET #", num2str(net_id)))
    disp('~~~~~~~~~~~~~~~~~~')
    [data_out,svnm] = SS_landscape_2D(N,base_params,param_names,S,Jmat,...
        p1,p2,pnum,p1min,p1max,p2min,p2max,net_id,fdr_nm);

    alldata_out(net_id) = {data_out};
    save(matname,'alldata_out')

    % Computer rests for 5 minutes after 20 iterations
    if c == 20
        pause(60*5)
        c = 0;
    end
    c = c + 1;
end


%% 3. Plot results of bifurcation / save workspace with plotting data

[p1_range,p2_range,all_valSSmap] = plot_2D_Bifirucation_v2(alldata_out,...
    p1min,p1max,p2min,p2max,pnum,p1,p2,param_names);

save(strcat('p',matname),'alldata_out','p1min','p1max','p2min','p2max',...
    'pnum','p1','p2','param_names','p1_range','p2_range','all_valSSmap')



