%% view_pat_info
%
% REQUIRED INPUTS
% pat_id
% run_nm (enter a workspace name if you want to avoid entering a bunch of
%   variables
%
% OPTIONAL INPUTS (only if run_nm = 'XXX.mat')
% run_mat
% TransMap
% SSnms
% data
% SSCounts
% TransMapMir
% SS_type
% NameMat
% PatSamNum
% TransCount

function [data4plotNaN, indivRunMat] = view_pat_info(pat_id,run_nm,run_mat,TransMap,SSnms,data,sp,ep,UID,SSCounts,TransMapMir,SS_type,NameMat,PatSamNum,TransCount,includeNaN)
    
    % for loading a workspace rather than input variables
    if nargin < 3
        load(run_nm,'TransMap','SSnms',...
            'data','SSCounts','TransMapMir','SS_type','NameMat','PatSamNum',...
            'TransCount','run_mat','sp','ep','UID')
        run_nm = extractBetween(run_nm,'_','.mat');
        includeNaN = false;
    elseif nargin > 2 && nargin < 16
        includeNaN = false; 
    end
    
    singleMap = squeeze(TransMap(pat_id,:,:));
    
    f1  = figure(1);
    subplot(4,2,[1 4])
    h = heatmap(singleMap);
    h.XDisplayLabels = SSnms(1:end-1);
    h.YDisplayLabels = SSnms(1:end-1);
    title(strcat(extractBefore(UID(sp(pat_id)),'_'),' - Transition Map'))
    ylabel('From/Start')
    xlabel('To/End')

    subplot(4,2,[7 8])
    bar(SSCounts(pat_id,:))
    xticks(1:7)
    xticklabels(SSnms(1:end-1))
    xtickangle(90)
    ylabel('Counts')
    title(strcat(extractBefore(UID(sp(pat_id)),'_'),' - SS Counts'))

    filenm = strcat(run_nm,'-',extractBefore(UID(sp(pat_id)),'_'),'-SS-Counts','.fig');
    saveas(f1,filenm)

    % Convert TransMap Mirror into

    singleMap = squeeze(TransMapMir(pat_id,:,:));
    Tidx = triu(ones(7,7)) == 1;

    TMap = singleMap;
    TMap(Tidx) = NaN;


    MirrorVectNaN = reshape(TMap,[],1);
    MirrorVect = MirrorVectNaN;
    MirrorVect(isnan(MirrorVectNaN)) = [];

    TransNmsNaN = reshape(NameMat,[],1);
    TransNms = TransNmsNaN;
    TransNms(isnan(MirrorVectNaN)) = [];

    % Rank the Most Frequenct (non-diagonal)
    [vals, idx] = maxk(MirrorVect,10);

    disp(['################################'])

    array2table([vals vals/PatSamNum(pat_id)],'RowNames',TransNms(idx),...
        'VariableNames', {'Counts','Frequency'})
    disp(strcat(extractBefore(UID(sp(pat_id)),'_'),' - Transition Info'))
    disp(['Total Samples: ', num2str(PatSamNum(pat_id))])
    disp(['Number of Transitions: ' num2str(TransCount(pat_id)), ' (', ...
        num2str(round(TransCount(pat_id)/PatSamNum(pat_id),3)),')'])

    % Plot Transitions alongside clinical data
    stress = str2double(data.STRESS);
    menses = ~isnan(str2double(data.MENSTRUATION));
    PH = str2double(data.PH);
    VagSex = str2double(data.VAG_INT);
    CST = str2double(data.CST_HL);

    IndivStrNaN = stress(sp(pat_id):ep(pat_id));
    IndivMensNaN = menses(sp(pat_id):ep(pat_id));
    IndivPHNaN = PH(sp(pat_id):ep(pat_id));
    IndivVagSexNaN = VagSex(sp(pat_id):ep(pat_id));
    IndivMapNaN = SS_type(sp(pat_id):ep(pat_id));
    IndivCSTNaN = CST(sp(pat_id):ep(pat_id));

    nanidx = IndivMapNaN == 8;

    data4plotNaN = [IndivMapNaN IndivStrNaN IndivMensNaN IndivPHNaN IndivVagSexNaN IndivCSTNaN];
    
    if includeNaN
        data4plot = data4plotNaN;
        SScols = [145 25 25;
        201 201 201;
        19 235 213;
        224 170 150;
        80 14 201;
        114 157 194;
        247 205 64;
        1 1 1]./255;
    else
        data4plot = data4plotNaN;
        data4plot(nanidx,:) = [];
        SScols = [145 25 25;
        201 201 201;
        19 235 213;
        224 170 150;
        80 14 201;
        114 157 194;
        247 205 64]./255;
    end




    f2 = figure(2);
    
    % Plot colorblcok/bin trajectory
    subplot(7,1,[1 2])
    imagesc(data4plot(:,1)')
    ylabel('SS Type')
    colormap(SScols)
    caxis([1 size(SScols,1)])
    cm = colorbar('southoutside');
    cm.TickLabels = SSnms(1:size(SScols,1));
    hold on
    cs = find(data4plot(:,2) > 3);
    plot(cs,ones(size(cs)),'o','MarkerEdgeColor',[1 1 1])
    cs5 = find(data4plot(:,4) > 5);
    plot(cs5,0.75.*ones(size(cs5)),'*','MarkerEdgeColor',[1 1 1])
    cs4 = find(data4plot(:,6) == 4);
    plot(cs4,1.25.*ones(size(cs4)),'v','MarkerEdgeColor',[1 1 1])
    me = find(data4plot(:,3) == 1);
    plot(me,1.1.*ones(size(me)),'d','MarkerEdgeColor',[1 1 1])
    title(strcat(extractBefore(UID(sp(pat_id)),'_'),' - Clincal Map'))

    % Plot actual relative abundance trajectory
    subplot(7,1,3)
    sp_cols = [0.9290 0.6940 0.1250;
        0.5 0.5 0.5;
        0.3010 0.7450 0.9330];
    sp_names = {'BV','LI','oLB'};
    
    indivRunMat = run_mat(sp(pat_id):ep(pat_id),:);
    colororder(sp_cols)
    for i = 1:3
        y = indivRunMat(:,i);
        x = 1:length(y);
        plot(x(~isnan(y)),y(~isnan(y)),'.:','MarkerSize',20)
        hold on
    end
    legend(sp_names)
    xlabel('Time (Days)')
    ylabel('Relative Abundance')

    % Plot CST
    subplot(7,1,4)
    plot(0.5:size(data4plot,1),data4plot(:,6),'Color','k')
    xlim([1 size(data4plot,1)])
    ylim([1 5])
    ylabel('CST - v')

    % Plot Stress Levels
    subplot(7,1,5)
    plot(0.5:size(data4plot,1),data4plot(:,2),'Color','k')
    xlim([0.5 size(data4plot,1)])
    ylim([1 5])
    ylabel('Stress - O')

    % Plot Menses
    subplot(7,1,6)
    plot(0.5:size(data4plot,1),data4plot(:,3),'Color','k')
    xlim([0 size(data4plot,1)])
    ylim([0 1.2])
    ylabel('Menses - d')

    
    % Plot PH
    subplot(7,1,7)
    plot(1:size(data4plot,1),data4plot(:,4),'Color','k')
    xlim([0 size(data4plot,1)])
    ylim([min(data4plot(:,4)) max(data4plot(:,4))+1])
    ylabel('PH - *')
%     plot(0.5:size(data4plot,1),data4plot(:,2))
%     xlim([0 size(data4plot,1)])
%     ylim([0 1])
%     ylabel('Vaginal Sex')

    filenm = strcat(run_nm,'-',extractBefore(UID(sp(pat_id)),'_'),'-Clincal-Map','.fig');
    saveas(f2,filenm)
    
end