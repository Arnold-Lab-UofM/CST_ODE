%% generate_Global_ParamSets.m
% [LHSmat,param_names,mat,mat_names] = generate_Global_ParamSets(all_params,std_params,nsample,sclF,logUniform,PID)
%
% Generates LHS parameter sets given set specifications.
%
% REQUIRES:
%   * lhs_ode_unif_new.m (Kirschner Global Sensitivity Analysis Code)
%   * get_info_SS_3sp.m
%
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% Mar 12, 2022
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [LHSmat,param_names,mat,mat_names] = generate_Global_ParamSets(all_params,std_params,nsample,sclF,logUniform,PID)
%%
sp_names = {'BV','Li','oLB'};
    
    num_sp = length(sp_names);
    
    [nm_out1] = generate_parameter_names(sp_names);
    [nm_out2] = generate_coeff_labels('\alpha',sp_names);

    param_names = horzcat(nm_out1{1:length(sp_names)},nm_out2);
    

    disp('###### Generated Uniform Dist. LHS ######')
    LHSmat = genUniform(nsample,all_params,std_params,sclF,num_sp,logUniform);
    ff = strcat('Unif-',num2str(sclF),'std');

   %% 
    Type = 'ralpha';
    [S, Jmat] = symbolic_solns(num_sp,Type); % calculate symbolic solution
    
    StbleSS = cell(size(LHSmat,1),1);
    parfor i = 1:size(LHSmat,1)
        params = LHSmat(i,:);
        [StableStates,~,~,~] = calc_SS_stability(num_sp,params,S,Jmat,Type);
        StbleSS{i} = {StableStates};
        disp(['Run #', num2str(i)])
    end

    % Plot Result
    [~,mat, ~,mat_names,~,~] = get_SS_info_3sp(StbleSS,false);
    
    out_nm = strcat(strcat(date,'-Global_summary-',num2str(sclF),'x-',...
        PID,ff,'.fig'));
    % savefig(gcf,char(out_nm))
    
    save(char(strcat(extractBefore(out_nm,'.fig'),'.mat')),'LHSmat',...
        'param_names','mat','mat_names','StbleSS','Type','S','Jmat')

end

function s = genNorm(nsample,all_params,std_params,num_sp)
    base_params = mean(all_params);
    s = NaN(nsample,length(base_params));
    
    gr_idx = 1:num_sp;
    si_idx = num_sp+1:num_sp+1:(num_sp^2 + num_sp);

    for i = 1:length(base_params)
        xmean = base_params(i);
        xsd = std_params(i);
        s(:,i) = lhs_ode_norm_new(xmean, xsd, nsample);
    end
    
    % Remove negative growth rates
    s(s(:,gr_idx) < 0) = 1E-10;
    
    % Remove positive self-interaction terms
    s(s(:,si_idx) > 0) = -1E-10;
end

function [s] = genUniform(nsample,all_params,std_params,sclF,num_sp,logUniform)
    
    base_params = all_params;
    s = NaN(nsample,length(base_params));
    
    gr_idx = 1:num_sp;
    si_idx = num_sp+1:num_sp+1:(num_sp^2 + num_sp);

    param_array = zeros(length(base_params),1);
    
    if sclF == 0
        param_array(:,1) = min(all_params);
        param_array(:,2) = max(all_params);
        
    else
        param_array(:,1) = base_params - sclF.*std_params;
        param_array(:,2) = base_params + sclF.*std_params;
        
    end

    
    param_array = sort(param_array,2);

    for i = 1:length(base_params)
        xmin = param_array(i,1);
        xmax = param_array(i,2);
        s(:,i) = lhs_ode_unif_new(xmin, xmax, nsample, logUniform);
    end
end