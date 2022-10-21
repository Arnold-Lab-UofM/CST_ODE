%% get_suggested_ics.m
% Numerically determines initial conditions needed to approach each
% possible steady-state in a multiple stable system.
%
% FULL ENTRY:
% [icTable, icmat, SS_type, SS_typenms] = get_suggested_ics(StableStates,...
%    base_params,nruns,tend,ODEoptions)
%
% OPTIONAL ENTRY:
% [icTable, icmat, SS_type, SS_typenms] = get_suggested_ics(StableStates,base_params)
% Here, default values for tend, ODEoptions and nruns are used.
%
% REQUIRED INPUTS:
% StableStates - m x n matrix of possible steady states for the parameter
%   set (m = number of steady states, n = number of species)
%
% base_params - your parameter set for the model
% 
% OPTIONAL INPUTS:
% nruns - number of ICs that are tried
% tend - how long the ODE solver runs
% options - the options variable for the ODE solver.
%
% OUTPUTS:
% icTable - Table that summarizes suggested ICs for each steady-state
%    possible
% icmat - full matrix of ICs that were used (nruns x num_sp matrix)
% SS_type - array that describes the SS_type for each IC in icmat
% SS_typenms - array that describes the name for each SS in SS_type arrary


function [icTable, icmat, SS_type, SS_typenms] = get_suggested_ics(StableStates,...
    base_params,nruns,tend,ODEoptions)

    % If user does not input nruns, tend, ODEoptions, use these defaults
    if nargin < 3
        tend = 500; % run for 500 hrs
        tspan = [0 tend];
        options = []; % default ODE solver options
        nruns = 50; % try 50 ICs
    else
        tspan = [0 tend];
        options = ODEoptions;
    end

    xmin = 0.001;
    xmax = 2*max(StableStates,[],'all');
    logUnif = true;
    icmat = [lhs_ode_unif_new(xmin, xmax, nruns, logUnif),...
        lhs_ode_unif_new(xmin, xmax, nruns, logUnif),...
        lhs_ode_unif_new(xmin, xmax, nruns, logUnif)];
    

    SS_type = ones(size(icmat,1),1);
    SS_typenms = cell(size(icmat,1),1);

    for i = 1:size(icmat,1)
        y0 = icmat(i,:);
        [t, y] = ode45(@lhs_ode_gLV,tspan,y0,options,base_params);
        yend = y(end,:);

        dif = sum((StableStates - repmat(yend,size(StableStates,1),1)).^2,2);
        [val,id(i)] = min(dif);

        relAbund = yend/sum(yend);

        [SS_type(i), SS_typenms(i), ~] = bin_to_SS(relAbund,false);
    end

    nms = {};
    tmp = unique(SS_type);
    
    disp('########## TRENDS IN IC DEPENDENCE ##########')
    if length(tmp) > size(StableStates,1)
        disp(['WARNING: Additional SS found, check length of simulation'])
    elseif length(tmp) < size(StableStates,1)
        disp(['WARNING: Less SS than expected found, run additional runs'])
    end

    ref_y0 = NaN(length(tmp),size(StableStates,2));
    for i = 1:length(tmp)
        idt = find(SS_type == tmp(i));
        ref_y0(i,:) = icmat(idt(1),:);

        nms{end + 1} = SS_typenms{idt(1)};

    end

    icTable = array2table(ref_y0,'RowNames',nms,'VariableNames',{'BV','LI','oLB'});
    
    
    disp(icTable)
%     Mdl = fitctree(icmat./sum(icmat,2),SS_typenms,'PredictorNames',{'BV','LI','oLB'});
%     view(Mdl,'Mode','Graph')
end