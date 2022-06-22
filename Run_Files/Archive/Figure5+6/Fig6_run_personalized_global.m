%% Fig6_run_personalized_global.m
%
% REQUIRES:
%   * LHS_best_params.mat
%   * generate_Global_ParamSets.m
%       * lhs_ode_unif_new.m (Kirschner Global Sensitivity Analysis Code)
%       * get_info_SS_3sp.m
%
% OUTPUT:
%   * Personalized steady-state space 
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% Mar 12, 2022
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% 1. Load input workspace (generated in Fig5_LHS_Clinical_Calibration.m)
fdr_loc = '../workspaces/'; 
load(strcat(fdr_loc,'LHS_best_params.mat'))
PID = PID_list;
%% 2. Run code to explore personalized parameter spaces using LHS

nsample = 1000; % Number of stratifications 
sclF = 0.50; % Bounds of LHS parameter ranges (base_parameter +/- sclF*base_parameter)
logUniform = false; % Whether ample from uniform distribution

id_vals = find(~isnan(all_fitmetG) == 1);

clear all_mat
for idx = 1:length(id_vals)
    all_params = best_paramset(id_vals(idx),:);
    std_params = all_params;
    disp('~~~~~~~~~~~~')
    disp(strcat('RUN FOR:', PID(idx)))
    
    [LHSmat,param_names,mat,...
        mat_names] = generate_Global_ParamSets(all_params,...
        std_params,nsample,sclF,false,PID(idx));
    close 
    
    all_mat(idx,:) = sum(mat)./sum(sum(mat));
end

%% 3. Visualize Result
h = heatmap(all_mat);
h.XDisplayLabels = mat_names;
h.YDisplayLabels = PID(id_vals);
title(strcat('Global-',num2str(nsample),'-samples-',num2str(sclF),'std-Unif'))
