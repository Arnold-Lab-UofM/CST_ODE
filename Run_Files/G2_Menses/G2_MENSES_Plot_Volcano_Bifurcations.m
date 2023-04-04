%% G2_MENSES_Plot_Compile_TimeSeries.m
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Goal: Plot comparisons of mono- and multi-stable equilibrium behaviors
% for populations simulated in the menses perturbation analysis.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Run after G2_MENSES_Run_TimeSeries_Analyis.m and
% G2_MENSES_Run_2D_Bifurcation.m
%   - Input: need to indicate location of the saved outputs for both of the
%   listed filles
%
% REQUIRES: plot_Volcano.m, plot_only_2D_Bifurcation.m
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This code plots the results for the volcano and bifurcation analyses
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 0. Set up locations of files to read in
clear;clc;
fd_id = 2; % 1 = oLB analysis, 2 = Li analysis
class_list = {{'1SS oLB', '2SS nAB or oLB'};
    {'1SS Li', '2SS nAB or Li'}};

% Enter locations of workspaces 
fdr_loc = '../G1_EB/';
load(strcat(fdr_loc,'SSConfig-Analysis-HMP.mat'),'S','Jmat','param_names')
fdr_names = strcat('MensesTimeSeries/',{'oLB-States-HMP_-446lhs-day7-29-Mar-2023';
    'Li-States-HMP_-309lhs-day7-29-Mar-2023'});
bif_ws = {strcat('MensesVolcanoBifurcation/',{'2D_BIF_1SS-oLB_HMP.mat','2D_BIF_2SS-oLB_HMP.mat'});
    strcat('MensesVolcanoBifurcation/',{'2D_BIF_1SS-Li_HMP.mat','2D_BIF_2SS-Li_HMP.mat'})};
ws_nm = '29-Mar-2023-4param-mod-0  0  0  0-for-7d-run.mat';
loc_name = strcat(fdr_names{fd_id},'/',ws_nm);

load(loc_name,'sel_nets','all_nm_CST')
N = 3; % number of species
if ~exist('all_nm_CST','var')
    all_nm_CST = [];
    for idx = 1:size(sel_nets,1)
        [run_mat,~,~,~] = calc_SS_stability(N,sel_nets(idx,:),S,Jmat);
        [nmf] = get_VALENCIA_class(run_mat); 
        all_nm_CST = [all_nm_CST; nmf];
    end
    save(loc_name,'all_nm_CST',"-append")
end

%% 1. PLOT VOLCANO
SS_names_CST = unique(all_nm_CST);

sel_idx = {};
for i = 1:length(SS_names_CST)
    tmp = all_nm_CST == SS_names_CST(i); % all_nm has been updated to all_nm_CST in some versions of code (Number of LHS samples x 1 string array)
    sel_idx{end+1} = find(tmp);
end
sel_nets1 = sel_nets(sel_idx{1},:);
sel_nets2 = sel_nets(sel_idx{2},:);

alpha = 0.01;
offset = 0.1;

f = figure;
plot_Volcano(sel_nets1,sel_nets2,alpha,offset,param_names,class_list{fd_id})

ftit = strcat(extractAfter(class_list{fd_id}{1}," "),"-menses_volcano.fig");
savefig(f,ftit)
close;
%% 2. PLOT BIFURCATIONS

bif_list = bif_ws{fd_id};
for k = 1:length(bif_list)
    ftit = strcat(class_list{fd_id}{k},"-menses_bifurcation.fig");
    load(strcat(bif_list{k}))
    f = figure;
    plot_only_2D_Bifurcation(all_valSSmap,p1_range,p2_range,p1,p2,param_names)
    savefig(f,ftit)
    close;
end

%% 3. Compile Manuscript Figure

fig_folder = 'MensesVolcanoBifurcation/'; % location of figures
listing = dir(fig_folder);
file_nms = {listing.name};
fig_nms = file_nms(contains(file_nms,'.fig'));

sp_id = contains(fig_nms,'oLB'); % pulls evaluation point
simType = contains(fig_nms,'volcano'); % keeps track of what type off simulation

volc_plots = [fig_nms(simType & sp_id), fig_nms(simType & ~sp_id)];
bif_plots = [fig_nms(~simType & sp_id), fig_nms(~simType & ~sp_id)];

fig_list = strcat(fig_folder,[volc_plots,bif_plots]);
idxs = {1,4,2,3,5,6};
nrows = 2;
ncols = 3;

finalfigure = combine_plots_to_subplot(fig_list,idxs,nrows,ncols);
savefig(finalfigure,strcat('compiled-',extractBefore(fig_folder,'/'),'-results.fig'))

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

