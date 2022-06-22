%% run_LHS_analysis_template.m

%% 1. Loads required workspace and sets up base variables
clear;clc;
fdr_loc = '../workspaces/'; 
load(strcat(fdr_loc,'SSConfig-Analysis-Model_LHS_10x.mat'), 'LHSmat','mat','StbleSS','S',...
    'Jmat','Type','SS_namesv','sp_names','param_names','all_nm')

% Prompt for SS type:
[indx,tf] = listdlg('PromptString',{'Pick Steady State Type:'},...
    'ListString',SS_namesv);
ss_id = indx;
[sel_idx] = find(all_nm == SS_namesv(indx));
sel_nets = LHSmat(sel_idx,:); %final parameter sets to use

%% 2. RUN SIMULATION  
pidx = [1]; % index of parameter of interest
% Lin = [1 7 14 21 60]; % duration
% var_vals = [-3 -2 -1 -0.75 -0.5 0]; % dose

Lin = [1 7];
var_vals = [-3 -2];


indx = 1; sp_idx = indx; % idx of dominant species
dom_ab = 0.7; min_ab = (1 - dom_ab)/2; ybase = repmat(min_ab,[1,3]);
ybase(sp_idx) = dom_ab;
plotTraj = false; % Generate individual plots (can take a long time to run)
time_post = 50; % how far after alteration to analyze
perChange = false; % relative to original parameter
ss_type = strrep(SS_namesv{ss_id},'/','');

for v = 1:length(var_vals)
    p1 = [var_vals(v)]; % Change in growth rate 1
    vectorCell = {p1}; % Get combinations
    for i = 1:length(Lin)
        L = Lin(i);
        sp_p = 7; % Start time of modification
        ep_p = sp_p + L; % End time of modification
        [arm] = simulate_LHS_time_response(ybase,ss_type,StbleSS,sel_nets,S,Jmat,Type,perChange,sp_p,ep_p,time_post,plotTraj,vectorCell,pidx);
    end
end


%% 3. Plot Results
dose = var_vals;
duration = Lin;
ev_point = 0; % Time point for evaluation (post regimen)
fldrnm = '2SS: [NO] CST-IV or [Li] CST-III_1D_k_grow-BV_22-Jun-2022';
plot_dosedur(fldrnm,dose,duration,ev_point)
