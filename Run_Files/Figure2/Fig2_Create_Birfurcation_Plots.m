%% Fig2_Create_Bifurcation_Plots.m
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% REQUIRED INPUTS:
%   * Model_LHS_5000_w_variables.mat
%   * ss_id: steady-state configuration (1SS - Li, 2SS - Li or BV...)
%   * net_id: selected representative parameter set
%
% REQUIRED FUNCTIONS
%   * SS_landscape.m 
%       * calc_SS_stability
%       * get_SS_info_3sp.m
%   * plot_colorblock_landscape.m
%
% NOTE: GUI prompts will guide user for entries.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% SETS USED IN THE MANUSCRIPT:
% (1) 1SS - oLB: ss_id = 7, net_id = 10
% (2) 2SS - oLB or BV: ss_id = 2, net_id = 16
% (3) 1SS - Li: ss_id = 6, net_id = 13
% (4) 2SS - Li and BV: ss_id = 3, net_id = 26
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% Jan 12, 2021
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%% 1. GENERATE OR LOAD REQUIRED DATA

% ####### MODIFY HERE ########
ss_id = 6; % steady-state configuration
net_id = 13; % which set to pull as representative
% ####### MODIFY HERE ########

fdr_loc = '../workspaces/';
load(strcat(fdr_loc,'Model_LHS_5000_w_variables.mat'),'LHSmat',...
    'mat', 'num_sp', 'S','Jmat','Type','colors','mat_names', ...
    'param_names')

sel_nets = LHSmat(mat(:,ss_id),:);
base_params = sel_nets(net_id,:);


%% 2. CALL 2D BIFURCATION FUNCTION
disp('~~~~~~~~~~~~~~~~~~')
disp(strcat("RUNNING: ", mat_names{ss_id}, " SET #", num2str(net_id)))
disp('~~~~~~~~~~~~~~~~~~')

clear SS_map data_out sum_table
[SS_map,data_out,sum_table,svnm] = SS_landscape(num_sp,base_params,...
    param_names,S,Jmat,Type,colors,ss_id);

%% 3. LOAD GENERATED WORKSPACE AND CREATE PLOT
plot_colorblock_landscape(svnm)

%% NOTE: The Model_LHS_5000_w_variables.mat just has the following data:
%
%   * sp_cols: colors for each microbial species in plots
%   * sp_names: names of species in manuscript
%   * param_names: parameter names for the model
%   * Steady-state configuration names and information:
%       * poss_SS: Counts of each of the 8 SS in model
%       * poss_SSnames: names for each SS in model
%       * mat_names: names of SS configurations in model
%       * mat: matrix of 0's and ones to denote which configuration per
%           parameter set
%       * mat_num: index of mat_names for corresponding steady state
%           configuration
%   * colors: colors for SS configuration landscape
%
%% CODE TO GENERATE FROM SCRATCH
% sp_cols = [0.9290 0.6940 0.1250;
%     0.5 0.5 0.5;
%     0.3010 0.7450 0.9330];
% 
% [poss_SS, mat, poss_SSnames, mat_names,mat_num,mat_code] = get_SS_info_3sp(StbleSS,false);
% 
% sp_names = {'BV','LI','oLB'};
% [nm_out1] = generate_parameter_names(sp_names);
% [nm_out2] = generate_coeff_labels('\alpha',sp_names);
% 
% param_names = horzcat(nm_out1{1:length(sp_names)},nm_out2);
% num_sp = 3;
%
% colors = [255 242 204;
%     193 147 171;
%     108 67 197;
%     180 198 231;
%     139 23 24;
%     194 194 194;
%     21 225 204;
%     168 99 255;
%     174 164 255;
%     255 183 255;
%     255 244 111;
%     231 0 255;
%     0 2 255;
%     224 170 150;
%     80 15 201;
%     114 157 195;
%     169 208 142;
%     84 130 53;
%     55 86 35;
%     191 143 0;
%     0 0 0]./255;

