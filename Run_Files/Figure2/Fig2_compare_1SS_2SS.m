%% Run File to Generate Figures in Figure 2
%
% GOAL: Determine parameters that most diffentiate two CST equilibirum
% types
%
% INPUT: 
%   - SSConfig-Analysis-Model workspace
%   - Select two CST Equilibrium Types to Compare
%
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 1. Load Data
clear;clc;
load('Figure1/SSConfig-Analysis-Model_LHS.mat')

%% 2. Pull Data (Choose 2 States to Compare)
[indx,~] = listdlg('ListString',SS_names_CST);


%% 3. Create volcano plot

classes = {SS_names_CST{indx(1)}, SS_names_CST{indx(2)}}; % labels for comparisons
sel_nets1 = LHSmat(all_nm == SS_names_CST{indx(1)},:);
sel_nets2 = LHSmat(all_nm == SS_names_CST{indx(2)},:);

alpha = 0.05; % significance threshold
offset = 0.05; % text offset

[PrismFormat] = plot_Volcano(sel_nets1,sel_nets2,alpha,offset,param_names,classes)