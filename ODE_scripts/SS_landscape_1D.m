%% [data_out,svnm] = SS_landscape_1D(num_sp,base_params,param_names,S,Jmat,...
%    pidx,pnum,pmin,pmax,net_id,fdr_nm)
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% INPUT:
%   * num_sp: number of species (3)
%   * base_params: model parameters to start the simulations
%   * param_names: names of the model parameters
%   * S: symbolic equations
%   * Jmat: Jacobian of ODE model
%   * bid: name of simulation
%   * pidx: parameter index
%   * pnum: number of steps
%   * pmin: minimum value in range
%   * pmax: maximum value in range
%   * net_id: reference index from LHS sets
%   * fdr_nm: where to save the files
%
% OUTPUT:
%   * data_out: predicted abundances for the steady states
%   * svnm: saved workspace name
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% v1: Jan 12, 2020
% v2: Jan 21, 2023 (Update to remove extra inputs, streamline for true 1
%   parameter global bifurcation)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [data_out,svnm] = SS_landscape_1D(num_sp,base_params,param_names,S,Jmat,...
    pidx,pnum,pmin,pmax,net_id,fdr_nm)
    
    p_range = linspace(pmin,pmax,pnum);
    num_iter = length(p_range);

    data_out = cell(num_iter,1);
    parfor i = 1:num_iter
        params = base_params;
        params(pidx) = base_params(pidx) + p_range(i);
        [StableStates,~,~,~] = calc_SS_stability(num_sp,params,S,Jmat);
        data_out{i} = {StableStates}; 
    end

    colors = [147	149	152;
        77	190	236;
        175	30	0;
        107	68	197;
        155	168	253;
        38	38	38;
        237	181	211;
        255	242	204;
        48	84	150;
        99	0	0;
        255	255	255]./255;

    % saves workspace
    op_nms = strrep(param_names,'\','');
    svnm = string(strcat(fdr_nm,'/N',num2str(net_id),'-',date,'-1D-',[op_nms{pidx}],'.mat'));
    svnm = strrep(svnm,'>','-');
    save(svnm,...
        'p_range', 'base_params', 'pidx','param_names','pmin','pmax','pnum', ...
        'data_out','colors')
end
