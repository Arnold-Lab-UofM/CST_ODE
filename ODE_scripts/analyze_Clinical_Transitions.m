%% [multi_stable_pid,multi_stable_types,mono_stable_pid, mono_stable_types] = analyze_Clinical_Transitions(fn)
%
% Function: Loads clinical data w/ transition matrices and determines
% multi-stability or mono-stability
%
% INPUT: 2022-04-29_VALENCIA-16s-to-SS-type.mat
%       - TransMap
%       - PatSamNum
%       - UID
%       - SSnms
%
% OUTPUT: Frequencies for multi-stability and mono-stability
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [multi_stable_pid,multi_stable_types,mono_stable_pid,...
    mono_stable_types] = analyze_Clinical_Transitions(fn)
    %% 1. Load Clinical Data Converted to Valencia CSTTs
    load(fn)

    %%
    prompt = {'Enter min samples:','Multi-Stable Threshold:',...
        'Enter Mono-Stable Threshold:','Enter Unstable Threshold:'};
    dlgtitle = 'Input';
    dims = [1 35];
    definput = {'10','0.2','0.5','0.2'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);

    L = str2double(answer{1});
    th1 = str2double(answer{2}); % Threshold of percent time points in a set state
    th = str2double(answer{3}); % Threshold for individuals in set state
    th2 = str2double(answer{4}); % Alternative state threshold (separates from multi-stable)


    %% 2. Extract patients with XX+ samples

    % Analyze how many patients have less than L time points
    less_10 = PatSamNum < L;
    PID_list = unique(extractBefore(UID,'_'));
    disp(strcat("### Patients with fewer than", num2str(L), " timepoints ###"))
    disp(['Number: ', num2str(sum(~less_10))])

    %% 3. ANALYZE TRANSITIONS: MultiStable
    bis = zeros(length(PatSamNum),1);
    typs = strings(length(PatSamNum),1);
    clear allq
    c = 1;
    for v = 1:length(PatSamNum)
        x = squeeze(TransMap(v,:,:))/(PatSamNum(v)-1); 
        q = diag(x);
        spl = q >= th1;
        %if isempty(
        %typs(v) = string(strcat(SSnms{find(spl)}));
        if sum(spl) == 2
            bis(v) = 1;
            typs(v) = string(strcat(SSnms{find(spl)}));
            %typs(c) = string(strcat(SSnms{find(spl)}));
            c = c + 1;
        end  
        allq(:,v) = q;
    end

    ms_idx = bis == 1 & PatSamNum > L & sum(allq)' > 0.6;
    % PID_list(ms_idx)
    disp(strcat("The MultiStable are: ", num2str(sum(ms_idx))))

    tabulate(typs)

    multi_stable_pid = ms_idx;
    multi_stable_types = typs(ms_idx);

    %% 4. ANALYZE STABLE INDIVIDUALS

    bis = zeros(length(PatSamNum),1);
    max_sp = zeros(length(PatSamNum),3);
    allq = zeros(3,length(PatSamNum));

    c = 1;
    for v = 1:length(PatSamNum)
        x = squeeze(TransMap(v,:,:))/(PatSamNum(v)-1); 
        q = diag(x);
        spl = q > th;
        va = sum(q(~spl) > th2);
        if sum(spl) == 1 && va == 0
            bis(v) = 1;

        end
        max_sp(v,:) = spl';
        allq(:,v) = q;
    end

    ss_idx = bis == 1 & PatSamNum > L; % get indexes
    PID_list(ss_idx)
    disp(strcat("The Single Stable are: ", num2str(sum(ss_idx))))

    disp("Breakdown by Type: [N0], [LI] or [oLB]")
    disp(sum(max_sp(ss_idx,:)))

    mono_stable_pid = ss_idx;
    mono_stable_types = max_sp(ss_idx,:);
end

