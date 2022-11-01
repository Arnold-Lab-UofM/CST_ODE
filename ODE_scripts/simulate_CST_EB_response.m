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

function [dirName] = simulate_CST_EB_response(ybase,ss_type,StbleSS,sel_nets,S,Jmat,Type,perChange,sp_p,ep_p,time_post,plotTraj,vectorCell,pidx)
 
    
    %% 1) Parameter information
    [~, mat, ~, mat_names,~] = get_SS_info_3sp(StbleSS,true);
    close
    sp_names = {'NO','LI','oLB'};
    [nm_out1] = generate_parameter_names(sp_names);
    [nm_out2] = generate_coeff_labels('\alpha',sp_names);
    param_names = horzcat(nm_out1{1:length(sp_names)},nm_out2);

    
    % Print the number of associated parameter sets
    disp('/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\')
    disp(['[',ss_type, '] has ', num2str(size(sel_nets,1)), ...
        ' parameter sets'])
    disp('\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/')

    %% 2. CREATE INPUT VALUE MATRIX FOR MODIFICATIONS
    
    if nargin < 8
        % Prompts user to select parameters to modify
        [pidx,~] = listdlg('PromptString', {'Select Parameters to Modify: '},...
            'ListString',param_names); 
    end

    if nargin < 7
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

    [newValueMat] = generate_input_combos(vectorCell);


    % Current "Defaults" - may change
    if nargin < 2
        perChange = false; % type of parameter alteration (absolute or relative)
        sp_p = 5; % Start time
        ep_p = sp_p + 72; % End time
        plotTraj = true; % Generate individual plots (can take a long time to run)
    end

    returnNorm = true; % type of parameter change (temporary or permanent)
    plotRel = 3; % close plots during simulation
    run_nets = sel_nets; % Input networks 

    % Display run information
    disp(['Altering SS Config. ', ss_type ,' run info: '])
    disp([param_names(pidx)])
    disp(['Duration of Alteration: ', num2str(ep_p - sp_p), 'd'])
    disp(['No. Parameter Sets: ',num2str(size(run_nets,1))])
    disp(['No. Parameter Combos: ', num2str(size(newValueMat,1))])
    disp('Input parameter values:')
    
    % Asker User to Confirm Run Details Before Proceeding
    str = input('Confirm Run: (Y/N)?','s');

    %5) START OF CODE NOT NEEDED TO BE MODIFIED
    if upper(str) == 'Y'
        tic
        pnm = strcat(param_names{[pidx]});
        pnmD = regexprep(pnm,{'{','}','\','>'},{'','','',''});
        dirName = strcat(ss_type,'_',num2str(length(pidx)),'D_',pnmD,'_',date);
        if not(isfolder(dirName))
            mkdir(dirName)
        else
            disp('WARNING: directory already exist - please modify name')
            newNm = inputdlg({'Enter New Name:'});
            dirName = newNm{1};
            mkdir(dirName)
        end

       
        disp('/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\')
        disp(['BEGINNING RUN: files to be saved in directory - ', dirName])
        disp('\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/')

        % Pre-allocate variables for parallel computing 
        sum_vals = NaN(size(newValueMat,1),3);
        allrunFlag = NaN(size(newValueMat,1),size(run_nets,1));
        allmaxChangeFlag = NaN(size(newValueMat,1),size(run_nets,1));
        cf = NaN(size(newValueMat,1));
        all_run_mat = cell(size(newValueMat,1),size(run_nets,1),1);
        NL = size(run_nets,1);

        % LOOP THRU ALL VALUES IN THE NEW VALUE MATRIX
        parfor val_id = 1:size(newValueMat,1)
            new_val = newValueMat(val_id,:);
            runType = cell(NL,1);
            runFlag = NaN(NL,1);
            maxChange = cell(NL,1);
            NewSS = cell(NL,1);
            maxChangeFlag = NaN(NL,1);

            % LOOP THROUGH ALL INPUT PARAMETER SETS
            for net_id = 1:NL
                base_params = run_nets(net_id,:);

                % Determine base SS
                [OGSS,~,~] = calc_SS_stability(length(sp_names),base_params,S,Jmat,Type);

                y0 = ybase*sum(OGSS(1,:)); % Set initial conditions based on input

                % START PLOTTING ALL TRAJECTORIES FOR GIVEN SET
                [tplot,yplot,newP,~] = change_parameter_men(base_params,y0,sp_p,ep_p,time_post,...
                    pidx,new_val,returnNorm,plotRel,perChange);
                title(['#',num2str(val_id), ': ',num2str(new_val)])
                run_mat = yplot;
                all_run_mat{val_id,net_id} = [tplot run_mat];

                % Determine altered SS
                [NewSS{net_id},~,~] = calc_SS_stability(length(sp_names),newP,S,Jmat,Type);

                % Get SS Classifications to determine Temp or Perm changes
                [SS_type, SS_typenms, SS_map] = bin_to_SS(run_mat,false);
                if SS_type(1) == SS_type(end)
                    runType{net_id} = 'TEMP';
                    runFlag(net_id) = 0; % Temp == 0
                    [~,popid] = max(y0); % maximum input is the starting SS
                    [minTraj,idx] = min(run_mat(:,popid)); % find where starting SS is a min (maximum effect)
                    maxChange{net_id} = SS_typenms{idx}; % save state of max effect
                    maxChangeFlag(net_id) = SS_type(idx);
                else
                    runType{net_id} = 'PERM'; % Perm classification may depend on length of simulation
                    runFlag(net_id) = 1;
                end
                disp(['****** New Val #' , num2str(val_id), ' - Set #', num2str(net_id), ' ******'])

            end

            % COMPILE INFO ON ALL THE RUNS
            allrunFlag(val_id,:) = runFlag;
            allmaxChangeFlag(val_id,:) = maxChangeFlag;

            % SUMMARY STATISTICS
            numPERM = sum(runFlag);
            if find(ybase == max(ybase)) == 1 % If starting in BV find Health
                numBVTEMP = sum(maxChangeFlag == 2 | maxChangeFlag == 3); % Reached Health temporarily 
                tmp = 'Health';
            else 
                numBVTEMP = sum(maxChangeFlag == 1); % Otherwise find BV
                tmp = 'BV';
            end
            numINTTEMP = sum(maxChangeFlag == 4 | maxChangeFlag == 5 | maxChangeFlag == 6 | maxChangeFlag == 7); % Had temporary change

            % DISPLAY RESULTS
            disp(['######## RUN: ', num2str(val_id), ' - ', num2str(new_val), ' ########'])
            disp(["Results at time: ", num2str(time_post), "Days after perturbation"])
            disp([num2str(numPERM), ' of ', num2str(length(runFlag)),...
                '(',num2str(round(numPERM/length(runFlag)*100,2)),...
                '%) of networks had permanent switch at'])
            disp([num2str(numBVTEMP), ' of ', num2str(length(runFlag)),...
            '(',num2str(round(numBVTEMP/length(runFlag)*100,2)),...
            '%) of networks had temp',tmp, 'switch'])
            disp([num2str(numINTTEMP), ' of ', num2str(length(runFlag)),...
            '(',num2str(round(numINTTEMP/length(runFlag)*100,2)),...
            '%) of networks had temp intermediate switch'])

            % COMPILE RESULTS
            sum_vals(val_id,:) = [numPERM numBVTEMP numINTTEMP];
        end

        if find(ybase == max(ybase)) == 1 % If starting in BV find Health
            tmp = 'Health';
        else 
            tmp = 'BV';
        end

        % SUMMARY TABLE NORMALIZED TO NUMBER OF RUNS
        summaryTable = array2table([newValueMat sum_vals./size(run_nets,1)*100], 'VariableNames',...
            horzcat(param_names(pidx),{'%Perm',strcat('%Temp',tmp),'%TempInt'}));
        [val, idx] = sort(summaryTable.("%Perm"),'descend');

        % PLOT SUMMARY STATS
        subplot(1,2,1)
        h1 = heatmap(summaryTable{idx,1:size(newValueMat,2)});
        h1.XDisplayLabels = summaryTable.Properties.VariableNames(1:size(newValueMat,2));
        h1.YDisplayLabels = idx;
        title(strcat('Parameter Values', ss_type))
        ylabel('Run Number')

        subplot(1,2,2)
        h2 = heatmap(summaryTable{idx,size(newValueMat,2)+1:end});
        h2.XDisplayLabels = summaryTable.Properties.VariableNames(size(newValueMat,2)+1:end);
        h2.YDisplayLabels = idx;
        title(strcat("System Trends @ ", num2str(time_post), "Days post Perturbation"))
        ylabel('Run Number')

        f = gcf;
        saveas(f,strcat(dirName,'/',strcat(dirName,'-','LHS-figure.fig')))

        % SAVE WORKSPACE AND WRITE SUMMARY TABLE TO EXCEL
        ws_nm = strcat(dirName,'/',date,'-',num2str(length(pidx)),'param-mod-',num2str(ep_p - sp_p), ...
            'hr-run.mat');
        save(ws_nm,'summaryTable','run_nets','sp_p','ep_p', 'pidx', 'newValueMat',...
            'param_names','ybase','allrunFlag','allmaxChangeFlag','vectorCell',...
            'returnNorm','plotRel','perChange','all_run_mat','dirName',...
            'ws_nm','time_post')
        % writetable(summaryTable,strcat(extractBefore(ws_nm,'.mat'),'.xls'))

        % PLOT ALL TRAJECTORIES
        if plotTraj
            close all; 
            LHS_trace_visualize_new(all_run_mat,newValueMat,dirName,pidx,sp_p,ep_p,time_post)
        end

        disp('/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\')
        disp(['RUN COMPLETE: saved in directory - ', dirName])
        disp('\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/')
        toc
    else
        disp('/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\')
        disp('RUN TERMINATED: Try again. Please enter Y or N')
        disp('\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/')
    end

%%%%%%%%%%%%%%%%%%%%%%%% END OF MAIN CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
