%% [dirName] = simulate_CST_EB_response.m
%
% REQUIRES: LHS_trace_analysis.m, generate_input_combos.m, and the parallel
% computing toolbox.
%
% INPUT (REQUIRED): name of workspace with parameter set information
%
% INPUT (OPTIONAL GROUP 1):
% - perChange: true if you are modifying a parameter relative to it's
%        value, false (default) if you are setting to a set new value
% - sp_p: starting point of alteration (default 5h)
% - ep_p: end point of alteration (default sp_p + 72h)
% - plotTraj: whether to plot each run and save the figure (default, true)
%
% INPUT (OPTIONAL GROUP 2):
% - pidx: index of each parameter to be modified
% - vectorCell: cell array with parameter modification
%
% OUTPUT:
% Files will be saved to a directory named off of the (1) SS configuration
%   (2) parameters changed (3) current date. 

% Includes workspace, summary figure of the runs, and if plotTraj is true, 
% all the trajectories for each parameter combination in the NewValueMat.
%
% Christina Y. Lee (May 20, 2021)

function [dirName] = simulate_CST_EB_response(ybase,ss_type,sel_nets,S,Jmat,perChange,sp_p,ep_p,time_post,plotTraj,vectorCell,pidx)
    
    %Parameter information
    sp_names = {'NO','LI','oLB'};
    [nm_out1] = generate_parameter_names(sp_names);
    [nm_out2] = generate_coeff_labels('\alpha',sp_names);
    param_names = horzcat(nm_out1{1:length(sp_names)},nm_out2);

    % Print the number of associated parameter sets
    disp('/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\')
    disp(['[',ss_type, '] has ', num2str(size(sel_nets,1)), ...
        ' parameter sets'])
    disp('\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/')

    % 2. CREATE INPUT VALUE MATRIX FOR MODIFICATIONS (if not all are set)
    if nargin < 7
        % Prompts user to select parameters to modify
        [pidx,~] = listdlg('PromptString', {'Select Parameters to Modify: '},...
            'ListString',param_names); 
    end

    if nargin < 6
        % Prompts user to input parameter values to try
        answ = inputdlg(param_names(pidx),'Input New Vals', [1 50], ...
            repmat({'[0.1 0.2 0.3]'},1,length(pidx)));

        vectorCell = cell(1,length(pidx));
        for i = 1:length(pidx)
            vectorCell{i} = str2num(answ{i});
        end
    else
        answ = 'User Input values: check order matches parameter order';
    end

    % Generates all possible combinations of parameters
    if length(vectorCell) == 1
        newValueMat = vectorCell{:};

    else
        [newValueMat] = generate_input_combos(vectorCell);
    end

    returnNorm = true; % type of parameter change (temporary or permanent)
    plotRel = 3; % close plots during simulation


    % Display run information
    disp(['Altering SS Config. ', ss_type ,' run info: '])
    disp([param_names(pidx)])
    disp(['Duration of Alteration: ', num2str(ep_p - sp_p), 'd'])
    disp(['No. Parameter Sets: ',num2str(size(sel_nets,1))])
    disp(['No. Parameter Combos: ', num2str(size(newValueMat,1))])
    disp('Input parameter values:')
    

    %5) START OF CODE NOT NEEDED TO BE MODIFIED

    tic
    pnm = strcat(param_names{[pidx]});
    pnmD = regexprep(pnm,{'{','}','\','>'},{'','','',''});
%     dirName = strcat(ss_type,'_',num2str(length(pidx)),'D_',pnmD,'_',date);
    dirName = strcat(ss_type,'_','-',num2str(size(sel_nets,1)),'lhs-day',num2str(sp_p),'-',date);
    if not(isfolder(dirName))
        mkdir(dirName)
    else
        disp('WARNING: directory already exist - please modify name')

    end

   
    disp('/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\')
    disp(['BEGINNING RUN: files to be saved in directory - ', dirName])
    disp('\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/')

    % Pre-allocate variables for parallel computing 
    all_run_mat = cell(size(newValueMat,1),size(sel_nets,1),1);
    NL = size(sel_nets,1);

    % LOOP THRU ALL VALUES IN THE NEW VALUE MATRIX
    [~,mxidx] = max(ybase);
    parfor val_id = 1:size(newValueMat,1)
        new_val = newValueMat(val_id,:);
        NewSS = cell(NL,1);
        % LOOP THROUGH ALL INPUT PARAMETER SETS
        for net_id = 1:NL
            base_params = sel_nets(net_id,:);

            % Determine base SS
            [OGSS,~,~] = calc_SS_stability(length(sp_names),base_params,S,Jmat);
            [~,midxss] = max(OGSS,[],2);
            [~,id] = intersect(midxss,mxidx);
            if isempty(id)
                id = 1;
            end
            
            y0 = ybase*sum(OGSS(id,:)); % Set initial conditions based on input

            % START PLOTTING ALL TRAJECTORIES FOR GIVEN SET
            [tplot,yplot,newP,~] = change_parameter(base_params,y0,sp_p,ep_p,time_post,...
                pidx,new_val,returnNorm,plotRel,perChange);
            title(['#',num2str(val_id), ': ',num2str(new_val)])
            run_mat = yplot;
            all_run_mat{val_id,net_id} = [tplot run_mat];

            % Determine altered SS
            [NewSS{net_id},~,~] = calc_SS_stability(length(sp_names),newP,S,Jmat);
            disp(['****** New Val #' , num2str(val_id), ' - Set #', num2str(net_id), ' ******'])

        end
    end

    % SAVE WORKSPACE AND WRITE SUMMARY TABLE TO EXCEL
    ws_nm = strcat(dirName,'/',date,'-',num2str(length(pidx)),...
        'param-mod-',num2str(newValueMat(1,:)),'-for-',num2str(ep_p - sp_p), ...
        'd-run.mat');
    save(ws_nm,'sel_nets','sp_p','ep_p', 'pidx', 'newValueMat',...
        'param_names','ybase','vectorCell',...
        'returnNorm','plotRel','perChange','all_run_mat','dirName',...
        'ws_nm','time_post')

    % PLOT ALL TRAJECTORIES
    if plotTraj
        close all; 
        plot_CST_EB_response(all_run_mat,newValueMat,dirName,pidx,sp_p,ep_p,time_post)
    end

    disp('/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\')
    disp(['RUN COMPLETE: saved in directory - ', dirName])
    disp('\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/')
    toc


%%%%%%%%%%%%%%%%%%%%%%%% END OF MAIN CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
