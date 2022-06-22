%% [SS_map,data_out,sum_table] = SS_landscape_global_loop(num_sp,base_params,param_names,S,Jmat,Type,colors,bid,varargin)
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% INPUT:
%   * num_sp: number of species (3)
%   * base_params: model parameters to start the simulations
%   * param_names: names of the model parameters
%   * S: symbolic equations
%   * Jmat: Jacobian of ODE model
%   * Type: formulation of gLV
%   * colors: colormap for plot
%   * bid: name of simulation
%   * varagin
%
% OUTPUT:
%   * SS_map: steady states indicator for each combination
%   * data_out: predicted abundances for the steady states
%   * sum_table: summary of steady states reached
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% Jan 12, 2020
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [SS_map,data_out,sum_table,svnm] = SS_landscape_global_loop(num_sp,base_params,param_names,S,Jmat,Type,colors,bid,varargin)
    CST_state = {'1SS: [Li] CST-III';'1SS: [oLB] CST-I/II/V';'1SS: [NO] CST-IV';'2SS: [NO] CST-IV or [oLB] CST-I/II/V';'2SS: [Li] CST-III or [oLB] CST-I/II/V';'2SS: [Li] CST-III or [Li] CST-III';'2SS: [NO] CST-IV or [Li] CST-III';'3SS: [NO] CST-IV or [Li] CST-III or [oLB] CST-I/II/V';'2SS: [oLB] CST-I/II/V or [oLB] CST-I/II/V';'2SS: [NO] CST-IV or [NO] CST-IV'};
    if nargin <= 8
        [p1,~] = listdlg('PromptString',{'Select Parameter 1'},'ListString',...
            param_names, 'SelectionMode','multiple');
        [p2,~] = listdlg('PromptString',{'Select Parameter 2'},'ListString',...
            param_names, 'SelectionMode','multiple');

        prompt = {'Log (Y/N)?','min','max','number runs'};
        dlgtitle = strcat('Input Range -',param_names{p1});
        dims = [1 75];
        definput = {'N','-0.2','0.2','100'};
        a1 = inputdlg(prompt,dlgtitle,dims,definput);

        dlgtitle = strcat('Input Range -',param_names{p2});
        a2 = inputdlg(prompt,dlgtitle,dims,definput);

        if a1{1} == 'Y'
            p1_range = logspace(str2double(a1{2}),str2double(a1{3}),str2double(a1{4}));
        else
            p1_range = linspace(str2double(a1{2}),str2double(a1{3}),str2double(a1{4}));
        end

        if a1{2} == 'Y'
            p2_range = logspace(str2double(a2{2}),str2double(a2{3}),str2double(a2{4}));
        else
            p2_range = linspace(str2double(a2{2}),str2double(a2{3}),str2double(a2{4}));
        end
    else
        p1 = varargin{1};
        p2 = varargin{2};
        p1_range = linspace(varargin{4},varargin{5},varargin{3});
        p2_range = linspace(varargin{6},varargin{7},varargin{3});

        net_id = varargin{8};
    end
    
    fdr_nm = varargin{9};

tic
    N1 = length(p1_range);
    N2 = length(p2_range);
    disp([N1,N2])
    
    if num_sp == 4
        get_info = @get_SS_info_4sp;
    elseif num_sp == 3
        get_info = @get_SS_info_3sp;
    end

    data_out = cell(N1,N2);

    parfor i = 1:N1
        for j = 1:N2
            params = base_params;
            params(p1) = base_params(p1) + p1_range(i)*abs(base_params(p1));
            params(p2) = base_params(p2) + p2_range(j)*abs(base_params(p2));
            [StableStates,SSval,eigval,UnstableStates] = calc_SS_stability(num_sp,params,S,Jmat,Type);
            data_out{i,j} = {StableStates};
        end
    end

    % Analyzing data_out
    SS_map = NaN(size(data_out));
    for i = 1:N1
        input_mat = data_out(i,:)';
        [poss_SS, mat, poss_SSnames, mat_names,mat_num,mat_code] = get_info(input_mat,false);
        SS_map(i,:) = mat_num';
    end
    
    if sum(sum(isnan(SS_map)))
       SS_map(isnan(SS_map)) = 21;

    end

    obs_ss = unique(SS_map);
    chk_nan = obs_ss == 21;

    if sum(chk_nan) ~= 0
        obs_ss = obs_ss(~chk_nan);
        obs_ss_nms = horzcat(mat_names(obs_ss),{'NaN'})';
        sum_table = array2table([obs_ss;21],'RowNames',obs_ss_nms,'VariableNames',{'mat_index'});
    else
        obs_ss_nms = horzcat(mat_names(obs_ss))';
        sum_table = array2table([obs_ss],'RowNames',obs_ss_nms,'VariableNames',{'mat_index'});
    end

    % Plot legend

    % saves workspace
    op_nms = strrep(param_names,'\','');
    svnm = string(strcat(fdr_nm,'/N',num2str(net_id),'-',date,'-2D-',[op_nms{p1}],'-vs-',[op_nms{p2}],'.mat'));
    svnm = strrep(svnm,'>','-');
    save(svnm,...
        'p1_range', 'p2_range', 'base_params', 'p1', 'p2','param_names', ...
        'data_out','SS_map','colors','mat_names','sum_table','mat_code')
    
    toc
end
