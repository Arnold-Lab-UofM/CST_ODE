%% Run File to Generate Figures in Figure 2
%
% GOAL: Determine parameters that most diffentiate two CST equilibirum
% types
%
% INPUT: 
%   - SSConfig-Analysis-Model workspace
%   - Select two CST Equilibrium Types to Compare
%
% OUTPUT:
%   - Workspace with results
%   - Spreadsheet with most consistent selected parameters
%   - Figures for min MSE and 1SE models
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 1. Load Data
clear;clc;
load('../workspaces/SSConfig-Analysis-Model_LHS_10x.mat')

%% 2. Pull Data (Choose 2 States to Compare)
[indx,~] = listdlg('ListString',SS_namesv);
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

%% 3. RUN LASSO

num_sets = 500; % number of sets too run
xnames = param_names; % parameter names
num_sel = 12; % number of parameters to plot
alpha = 1; % LASSO = 1, Ridge = 0

% auto generate file names:
nm1 = strcat(SS_namesv(indx(1)),"_VS_",SS_namesv(indx(2)));
nm2 = regexprep(nm1,{':',' ','/'},{'-','',''});
fittl = strcat(nm2,'.mat');

[T_sorted_1SE, T_sorted_min, fc_count] = Parallel_Resample_BINOM(num_sets,xblock,yblock,xnames,num_sel,alpha,fittl);
