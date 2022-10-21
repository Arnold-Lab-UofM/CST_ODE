%% convert_HMP_to_SS_types(ws_name)
%
% Requires workspace that has entry "run_mat", HMP "data", patient names
% "UID", "sp", and "ep".
%
%

%%
function [fn] = convert_HMP_to_SS_types(ws_name,run_title,SSref)
    %% Load data
    load(ws_name,'run_mat','data','UID','sp','ep')
    
    %% Classify
    
    SSnms = {'Only BV', 'Only Li', 'Only Lc', ...
    'Bv and Li', 'Bv and Lc', 'Li and Lc', ...
    'Bv, Li and Lc','NaN'};

    if nargin < 3 % default "bins"
        SSref = [1.0 0.0 0.0;
            0.0 1.0 0.0;
            0.0 0.0 1.0
            0.5 0.5 0;
            0.5 0 0.5;
            0.0 0.5 0.5;
            0.33 0.33 0.33];
    end

    SS_type = ones(size(run_mat,1),1);
    SS_typenms = cell(size(run_mat,1),1);
    SS_map = zeros(size(run_mat,1),size(SSref,1));

    for i = 1:size(run_mat,1)
        tmp = run_mat(i,:);
        if sum(isnan(tmp)) > 0
            id = 8;
        else

            dif = sum((SSref - repmat(tmp,size(SSref,1),1)).^2,2);

            [val,id] = min(dif);
        end

        SS_type(i) = id;
        SS_typenms(i) = SSnms(id);
        SS_map(i,id) = 1;
    end

    %%
    T = [run_mat SS_type];

    Ts = array2table(T,'VariableNames',{'BV','LI','oLB','ClassNum'});
    Ts.('ClassNm') = SS_typenms;

    %% Prepare Naming Conventions
    NameMat = cell(7,7);

    for i = 1:size(NameMat)
        for j = 1:size(NameMat)
            NameMat{i,j} = strcat('[',SSnms{i},'] w [', SSnms{j},']');
        end
    end
    
    %% Loop Thru all Patients
    TransMap = zeros(length(ep),7,7);
    TransMapMir = zeros(length(ep),7,7);
    TransCount = zeros(length(ep),1);
    PatSamNum = zeros(length(ep),1);
    UniqSS = cell(length(ep),1);
    SSCounts = zeros(length(ep),7);

    for i = 1:length(ep)
        tmp_p = SS_type(sp(i):ep(i));
        rm_tmp = tmp_p;
        rm_tmp(tmp_p == 8) = [];
        PatSamNum(i) = length(rm_tmp);
        UniqSS(i) = {unique(rm_tmp)};
        SSCounts(i,rm_tmp(1)) = SSCounts(i,rm_tmp(1)) + 1;
        for j = 1:length(rm_tmp)-1
            id1 = rm_tmp(j);
            id2 = rm_tmp(j+1);
            SSCounts(i,id2) = SSCounts(i,id2) + 1;
            if id1 == id2
                TransMap(i,id1,id2) = TransMap(i,id1,id2) + 1;
            else
                TransMap(i,id1,id2) = TransMap(i,id1,id2) + 1;
                TransCount(i) = TransCount(i) + 1;
            end
        end
        TransMapMir(i,:,:) = squeeze(TransMap(i,:,:)) + squeeze(TransMap(i,:,:))';
    end

    %% Frequency of transitions
    FrqTrans = TransCount./PatSamNum;
    
    %% Save info into one variable
    fn = strcat(datestr(now,'yyyy-mm-dd'),'_',run_title,'-16s-to-SS-type.mat');
    save(fn, 'run_mat', 'FrqTrans', 'TransMapMir', 'SSCounts','UniqSS',...
        'PatSamNum','TransCount','TransMap', 'sp','ep', 'NameMat', ...
        'SSnms', 'UID','data','SS_type','SS_typenms','SSref','Ts')
end