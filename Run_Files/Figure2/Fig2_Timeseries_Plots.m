%% Fig2_Timeseries_Plots.m
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% REQUIREMENTS:
%   * change_parameter.m
%   * calc_SS_stability.m
%   * bin_to_SS.m
%   * Model_LHS_5000_w_variables.mat
%       % See Fig2_Creat_Bifurcation_Plots.m for explanation of this
%           workspace
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
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

%% 1. LOAD DATA AND PREP WORKSPACE

% ####### MODIFY HERE ########
ss_id = 6; % steady-state configuration
net_id = 13; % which set to pull as representative
% ####### MODIFY HERE ########
fdr_loc = '../workspaces/';
load(strcat(fdr_loc,'Model_LHS_5000_w_variables.mat'),'LHSmat',...
    'mat', 'num_sp', 'S','Jmat','Type','colors','mat_names', ...
    'param_names','sp_names')

sel_nets = LHSmat(mat(:,ss_id),:);
base_params = sel_nets(net_id,:);

%% 2. TEMPORARY CHANGE IN PARAMETER

% MODIFY THESE FOR RUN
sp_p = 24; % @24 hr change parameter
returnNorm = true; % revert parameter to original value
ep_p = 96; % @96 hr return parameter to normal
vals = [-0.013 -0.0125];
ybase = [0.1 2 0.1];
relPlot = true; % plot relative abundance
% END MODIFICATION

[pidx,tf] = listdlg('ListString',param_names); % Pick parameter to change
c = 1;
for new_val = vals
    y0 = [1 1 1].*ybase;
    
    figure(1)
    ax1 = subplot(1,length(vals),c);
    [tplot,yplot,newP,pidx] = change_parameter(base_params,y0,sp_p,ep_p,pidx,new_val,returnNorm,relPlot);
    [StableStates,~,~] = calc_SS_stability(length(sp_names),newP,S,Jmat,Type);
    b = gca; legend(b,'off');
    
    figure(2)
    ax2 = subplot(1,length(vals),c);
    run_mat = yplot./sum(yplot,2);
    [SS_type, SS_typenms, SS_map] = bin_to_SS(run_mat,true);
    title(['Change: ', param_names{pidx}, ' from ', num2str(round(base_params(pidx),3)), ...
        ' to ', num2str(newP(pidx))])
    hold on
        q = xline(sp_p,'--',{'Altered', param_names{pidx}});
    hold on
    
    if returnNorm 
        q2 = xline(ep_p,'--',{'Altered', param_names{pidx}});
    end
    
    xlim([sp_p-20 length(SS_type)])
    
    hold off

    c = c + 1;
end

