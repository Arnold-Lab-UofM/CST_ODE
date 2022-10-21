%% Fig2_compare_1SS_2SS.m
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% GOAL: Determine parameters that most diffentiate two CST equilibirum
% types. Comparisons in manuscript:
%   (1) 1SS: [Li] CST -III vs 2SS: [Li] CST -III or [NO] CST-IV
%   (2) 1SS: [oLB] CST -I/II/V vs 2SS: [oLB] CST -I/II/V or [NO] CST-IV
%   (3) 1SS: [NO] CST -IV vs 2SS: [Li] CST -III or [NO] CST-IV
%   (4) 1SS: [oLB] CST -I/II/V vs 2SS: [oLB] CST -I/II/V or [NO] CST-IV
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% INPUT: 
%   * SSConfig-Analysis-Model workspace
%   * Select two CST Equilibrium Types to Compare
%
% OUTPUT:
%   * Workspace with results
%   * Spreadsheet with most consistent selected parameters
%   * Figures for min MSE and 1SE models
%
% REQUIRES:
%   * Add General_scripts to path
%       * Parallel_Resample_BINOM.m
%           * LASSO_ELASTIC_BINOMIAL.m
%   * Parallel Computing Toolbox
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee (chyylee@umich.edu)
% University of Michigan
% Jun 22, 2022
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 1. Load Data
% All workspaces used in our analysis are in the workspaces folder

load('../workspaces/SSConfig-Analysis-Model_LHS_10x.mat')

%% 2. Pull Data (Choose 2 States to Compare)

% prompts user for EB behavior selection
[indx,~] = listdlg('ListString',SS_namesv); 

% Pulls and formats data for LASSO into xblock and yblock 
sel_id = NaN(size(all_nm,1),length(indx));
sel_ps = cell(length(indx),1);
xblock = [];
grp = [];
for i = 1:length(indx)
    sel_id(:,i) = all_nm == SS_namesv{indx(i)};
    tmp = LHSmat(sel_id(:,i)==1,:);
    sel_ps(i) = {tmp};   
    xblock = [xblock; tmp];
    grp = [grp;ones(size(tmp,1),1)*indx(i)];
end

for i = 1:length(indx)
    yblock(:,i) = grp == indx(i);
end

%% 3. Run Re-sampling LASSO
% This code runs LASSO on re-sampled data. The selected classes are
% re-sampled to minimize the effect of class size, and LASSO is run
% multiple times as results can vary depending on how the data is split.
% The results then give a degree of confidence on the impact of each
% variable in discriminating between the groups analyzed.

% Enter information for LASSO run
num_sets = 500; % number of sets too run
xnames = param_names; % parameter names
num_sel = 12; % number of parameters to plot
alpha = 1; % LASSO = 1, Ridge = 0

% Auto generates file names:
nm1 = strcat(SS_namesv(indx(1)),"_VS_",SS_namesv(indx(2)));
nm2 = regexprep(nm1,{':',' ','/'},{'-','',''});
fittl = strcat(nm2,'.mat');

delete(gcp('nocreate'))
% Call the LASSO code
[T_sorted_1SE, T_sorted_min, fc_count] = Parallel_Resample_BINOM(num_sets,...
    xblock,yblock,xnames,num_sel,alpha,fittl);
