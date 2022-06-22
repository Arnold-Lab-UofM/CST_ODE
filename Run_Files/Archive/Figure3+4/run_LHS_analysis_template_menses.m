%% run_LHS_analysis_template_menses.m
%
% Example code for the analysis in figures 3 and 4.
%
% Example of how to run the perturbation simulations, using menses as an
% example. The core code has the capibility to screen for many different
% combinations of parameter perturbations. Below, it is looking at how
% changes to 4 different parameters at varying degrees of parameter value
% alterations impacts community composition over time.
%
% REQUIRES:
%   * Model_LHS_5000_w_variables.mat
%   * simulate_LHS_response_men.m
%       * change_parameter_men.m
%       * LHS_trace_visualize_new.m
%   * check_comp (at end of scirpt)
%
% NOTE: Some GUI and command window prompts will appear. This is to help
% guide the user and ensure the correct inputs are selected.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% Jan 12, 2021
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% 1.) Loads required workspace and sets up base variables
fdr_loc = '../workspaces/'; 
load(strcat(fdr_loc,'Model_LHS_5000_w_variables.mat'), 'LHSmat','mat','StbleSS','S',...
    'Jmat','Type','mat_names','sp_names','param_names')

% Prompt for SS type:
[indx,tf] = listdlg('PromptString',{'Pick Steady State Type:'},...
    'ListString',mat_names);
ss_id = indx;

init_state = 'LB'; % Starting in LB dominant state
[sel_idx] = check_comp(ss_id,mat,StbleSS,init_state);

sel_nets = LHSmat(sel_idx,:); %final parameter sets to use

%% 2. Run the perturbation analysis (Below is an example for Menses)
%   To modify for ABX, select BV as the initial dominating species, change
%   perChange to false (now will subtract rather than fold change), change
%   pidx = 1 (kgrow_BV), p1 = 0.95 d^-1 (change in growth rate), and
%   vectorCell should now only have one entry {p1}, make sure to change
%   "init_state" above to 'BV'

% Prompts for dominate species (L. iners or oLB, depending on SS type)
[indx,tf] = listdlg('PromptString',{'Pick Initial Dominating Species:'},'ListString',sp_names);
sp_idx = indx; % idx of dominant species

% ###### MODIFY BELOW ######
L = 7; % Duration of chanve
sp_p = 7; % Start time of modification
ep_p = sp_p + L; % End time of modification
plotTraj = true; % Generate individual plots (can take a long time to run)
perChange = true; % true = fold change, false = subtract from original
pidx = [2,3,7,10]; % index of growth rate parameters for L. iners and oLB and interactions terms Li/oLB -> BV
p1 = [0.5 0.75]; % Fold-change in growth rate Li
p2 = p1; % Fold-change in growth rate oLB
p4 = [0 0.5 1]; % Fold-change in interaction of Li -> BV
p5 = p4; % Fold-change of interation of oLB -> BV
vectorCell = {p1,p2,p4,p5}; % Get combinations
time_post = 50; % Length of simulation after perturbation
dom_ab = 0.7; % setting up initial conditions
min_ab = (1 - dom_ab)/2; % setting up initial conditions
ybase = repmat(min_ab,[1,3]); % setting up initial conditions
ybase(sp_idx) = dom_ab; % setting up initial conditions
ss_type = mat_names{ss_id}; % run name (SS configuration)
% ###### MODIFY RUN ABOVE ######

% Call the function
simulate_LHS_response_men(ybase,ss_type,StbleSS,sel_nets,S,Jmat,Type,...
    perChange,sp_p,ep_p,time_post,plotTraj,vectorCell,pidx)

%% check_comp
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% check_comp checks that the input parameter set is associated with a
% certain composition at steady-state (BV or LB).
%
% INPUT:
%   * ss_id = index of desired SS configuration 
%   * mat = matrix to get rows indexes to pulled desired parameter sets
%   * StbleSS = cell array of analytically determine SS compositions
%   * init_state = whether you want to start with BV dominance of LB
%       dominance (enter 'BV' or 'LB')
%
% OUTPUT:
%   * sel_idx = final indexes after verifying composition
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [sel_idx] = check_comp(ss_id,mat,StbleSS,init_state)
    if contains(upper(init_state),'BV')
        th = 0.3;
        th_func = @(x,th) x < th;
    elseif contains(upper(init_state),'LB')
        th = 0.7;
        th_func = @(x,th) x > th;
    else
        disp('WARNING: Please enter BV or LB for init_state')
    end
    
    sel_idx = [];
    idxs = find(mat(:,ss_id));
    for j = 1:length(idxs)
        c = 0;
        tid = idxs(j);
        tmp = StbleSS{tid}{:};
        rtmp = tmp ./ sum(tmp,2);
        for k = 1:size(rtmp,1)
            x = sum(rtmp(k,[2,3]));
            if th_func(x,th) % checks dominance (70%)
                c = c + 1;
            end
        end   
        if c >= 1
            sel_idx = [sel_idx;tid];
        end        
    end
end
    