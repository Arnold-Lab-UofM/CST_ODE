%% change_abundance.m
% [tpre,ypre,tpost,ypost] = change_abundance(base_params,y0,sp_p,ep_p,sp_idx,sp_amt)
%
% Generates a simulation of changing the abundance of a species in the
% model (like adding a probiotic).
%
% NON-GUI RUN:
% [tpre,ypre,tpost,ypost] = change_abundance(base_params,y0,sp_p,ep_p,sp_idx,sp_amt)
%
% GUI PROMPTED RUN:
% [tpre,ypre,tpost,ypost] = change_abundance(base_params)
% 
% REQUIRED INPUT:
% base_params - parameter set for simulation
%
% OPTIONAL INPUT:
% y0 - initial conditions
% sp_p - start point for added species
% ep_p - end point of simulations
% 
% sp_idx - species to add (index)
% sp_amt - amount of species to add
%
% OUTPUT: trajectory info
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% Jan 12, 2021
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [tpre,ypre,tpost,ypost] = change_abundance(base_params,y0,sp_p,ep_p,sp_idx,sp_amt)

    %% Load base colors and names
    sp_cols = [0.9290 0.6940 0.1250;
    0.5 0.5 0.5;
    0.3010 0.7450 0.9330];
    sp_names = {'BV','LI','oLB'};
    num_sp = length(sp_names);
    
    %% If GUI is activated: (only enter base_params)
    if nargin < 2
        % call SS function
        Type = 'RALPHA';
        [S, Jmat] = symbolic_solns(num_sp,Type);
        [StableStates,~,~] = calc_SS_stability(num_sp,base_params,S,Jmat,Type);
        
        % PROMPT FOR IC SELECTION
        if size(StableStates,1) > 1 % if multi-stable will ask user to select base SS
            [icTable, ~, ~, ~] = get_suggested_ics(StableStates,...
                base_params,50,500,[]);

            [idx,tf] = listdlg('PromptString', {'Select base Steady-State'},...
                'ListString',icTable.Properties.RowNames);

            y0 = icTable{idx,:};
        else % otherwise the default y0 is all ones
            y0 = ones(1, num_sp);
        end
        
        % PROMPT FOR SPECIES
        [sp_idx,tf] = listdlg('PromptString', {'Choose sp. to add'}, ...
            'ListString',sp_names);
        
        % PROMPT FOR TIME POINTS AND AMOUNT
        prompt = {'Perturbation Time','Total Run Time:','Amount to Add:'};
        dlgtitle = 'Input';
        dims = [1 35];
        definput = {'24','144','1'};
        A = inputdlg(prompt,dlgtitle,dims,definput);
        
        sp_p = str2double(A{1});
        ep_p = str2double(A{2});
        sp_amt = str2double(A{3});
        
    end
    
    %% Run simulation
    tspan = [0 sp_p];
    options =[];
    [tpre, ypre] = ode45(@lhs_ode_gLV,tspan,y0,options,base_params);
    
    y0 = [ypre(end,:)];
    y0(sp_idx) = sp_amt;
    tspan = [sp_p ep_p];
    options =[];
    [tpost, ypost] = ode45(@lhs_ode_gLV,tspan,y0,options,base_params);
    
    colororder(sp_cols)
    plot([tpre;tpost],[ypre;ypost],'LineWidth',1.5)
    legend(sp_names)
    xlabel('Time (h)')
    ylabel('Abundance')
    set(gca,'fontsize',12)
    xlim([sp_p-20 ep_p])
    hold on
    q = xline(sp_p,'--',{'Added:', sp_names{sp_idx}});
    legend([sp_names,'Perturb'])
    title([sp_names{sp_idx}, ' ', num2str(sp_amt), ' added'])
end