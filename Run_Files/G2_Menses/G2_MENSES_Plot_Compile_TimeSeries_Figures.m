%% G2_MENSES_Plot_Compile_TimeSeries.m
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Goal: Visualize the results of the menses perturbation simulations and
% compare results to the clinical data of the HMP cohhort.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Run after G2_MENSES_Run_TimeSeries_Analyis
%
% REQUIRES: 
%   - superbar.m https://www.mathworks.com/matlabcentral/fileexchange/57499-superbar
%   - Folder name and location for G2_MENSES_Run_TimeSeries_Analysis
%       outputs
%   - Clinical HMP data with menses metadata (Mar16-2023-EB.mat)
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This code plots the results for all time series data.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% ENTER DETAILS OF SIMULATION
fd_id = 1; % Index of folder of interest (see fdr_names for order)
threshhold = 0.5; % threshold for species dominaance
sp_idx = 2;  % initial species dominant (1 - nAB, 2 - Li or oLB)
menses_mag_idx = 4; % what menses magnitude simulation (this is -2x/1x)

%% 1. INDICATE LOCATION OF MENSES RESULTS & NAMES OF FILES
output_folder = 'MensesTimeSeries/';
fdr_names = {'oLB-States-HMP_-446lhs-day7-29-Mar-2023';
    'Li-States-HMP_-309lhs-day7-29-Mar-2023'};

ws_nm = '29-Mar-2023-4param-mod-0  0  0  0-for-7d-run.mat';
loc_name = strcat(output_folder,fdr_names{fd_id},'/',ws_nm);

%% 2. Simulatioon Time Series Data
th = threshhold; % threshold for composition to switch
xlimit = 10; % how much time to plot post-menses

f1 = figure(1);
figtit = strcat(extractBefore(loc_name,'_'),'-Simulation-th',...
    strrep(num2str(th),'.',''),'-mag',strrep(num2str(menses_mag_idx)," ",""),date,'.fig');
[Counts,total_runs] = plot_multi_panel_menses(loc_name,th,sp_idx,menses_mag_idx,xlimit); % plotting multi panel function
SimulationCounts = [Counts(:,1) total_runs - Counts(:,1)];
savefig(f1,figtit)

%% 3. Clinical Time Series Data
% Load HMP (clinical) data
load('../workspaces/Mar16-2023-EB.mat')

clin_pairs = [3, 6; % 1SS oLB w/ 2SS oLB or nAB
    2,4]; % 1SS Li w/ 2SS Li or nAB

inp_mat = indexAssignedEB(allAssigned == clin_pairs(fd_id,1) | allAssigned == clin_pairs(fd_id,2));

dur_req = 3; % minimum menses length
tp = 5; % time pre/post to analyze menses

[alldurRM, allpostRM, allpreRM,allmaxdur] = pull_menses_info(inp_mat,...
    data,sp,ep,run_mat,dur_req,tp);

% Plot Data
ci_th = 0.95; % percent confidence internval
f2 = figure(2);
figtit = strcat(extractBefore(loc_name,'_'),'-Clinical-th',...
    strrep(num2str(th),'.',''),'.fig');

switch_ids = plot_clinical_menses(allpreRM,allpostRM,alldurRM,allmaxdur,tp,th,ci_th);
savefig(f2,figtit)
ClinicalCounts = [sum(switch_ids == 1), sum(switch_ids == 2)];

%% 4. Compare Simulation to Clinical Data

f3 = figure(3);
P = plot_Simulation_vs_Clinical(ClinicalCounts,SimulationCounts);
figtit = strcat(extractBefore(loc_name,'_'),'-Clinical_vs_Simulation-th',...
    strrep(num2str(th),'.',''),'-mag',strrep(num2str(menses_mag_idx)," ",""),'.fig');

savefig(f3,figtit)

%% 5. Compile All Plots For Manuscript Figure

if length(menses_mag_idx) == 1 % main text
    nrows = 2;
    ncols = 4;
    fig_list = [f1,f2,f3];
    idx = {[1 2 3 4],[5 6 7],8};
    sz = [1,1,8,3];
else
    nrows = length(menses_mag_idx)+1; % supplement
    ncols = 4;
    fig_list = [f1,f2,f3];
    idx = {[1:16],[17:19],20};
    sz = [1,1,8,8];
end

FF = combine_plots_to_subplot(fig_list,idx,nrows,ncols);

set(FF,'units','inches','position',sz,'renderer','painters')

figtit = strcat(extractBefore(loc_name,'_'),'-Compiled-th',...
    strrep(num2str(th),'.',''),'-mag',strrep(num2str(menses_mag_idx)," ",""),'.fig');
savefig(FF,figtit)

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Functions Specific to this File
%   * plot_Simulation_vs_Clinical
%   * pull_menses_info
%   * plot_clinical_menses
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function P = plot_Simulation_vs_Clinical(ClinicalCounts,SimulationCounts)
% Inputs: 
%   * ClinicalCounts [Counts for switch at day 0, Clinical Data]
%   * Simulation Counts [Counts for swich at day 0, Simulation]
%
% Outputs:
%   * Symmetric marix of P-values for all pairwise comparisons

    if SimulationCounts(1,1) == 0 % remove null runs
        CountList = [ClinicalCounts;
            SimulationCounts(2:end,:)];
    else
       CountList = [ClinicalCounts;
            SimulationCounts];
    end
    
    numComparisons = size(CountList,1);
    P = NaN(numComparisons,numComparisons);
    for i = 1:numComparisons
        for j = 1:numComparisons
            observed = [CountList(i,:)', CountList(j,:)'];
            [p, chi2stat,df] = chigof(observed);
            P(i,j) = p;
        end
    end
    
    
    Y = [CountList(:,1)./sum(CountList,2)*100];
    E = zeros(size(Y));
    
    superbar(Y, 'E', E, 'P', round(P,4),'BarFaceColor','b','BarEdgeColor','none')
    for i = 1:size(Y,1) % Label for percent nAB response post
        text(i,Y(i,1),num2str(Y(i,1),'%.1f %%'),'VerticalAlignment','top', 'HorizontalAlignment', 'center',...
            'fontname','Arial','Color','w')
    end
    ax = gca;
    ax.LineWidth = 1;
    ax.XColor = 'k'; % Red
    ax.YColor = 'k'; % Blue
    ax.FontSize = 7;
    ylabel('Percent Samples')
    xticks(1:numComparisons)
    xticklabels({'Clinical','Model1','Model2','Model3'})
end

%%
function [alldurRM, allpostRM, allpreRM,allmaxdur] = pull_menses_info(inp_mat,data,sp,ep,run_mat,dur_req,tp)
% Inputs: 
%   * inp_mat: indexes of clinical data to anallyze
%   * data: table array of raw data from HMP study
%   * sp: indexes of start point for each patient
%   * ep: indexes of the end points of each patient
%   * run_mat: relative abundances converted to nAB, Li and oLB
%   * dur_req: duration requirment for menses (days)
%   * tp: time period to assess before and after menses
%
% Outputs:
%   * alldurRM: run_mat during menses
%   * allpostRM: run_mat post menses
%   * allpreRM: run_mat pre menses
%   * allmaxdur: maximum values of nAB during menses

    c = 1;
    for k = 1:length(inp_mat)
        pat_id = inp_mat(k);
        menses = ~isnan(str2double(data.MENSTRUATION));
    
        indivMenses = menses(sp(pat_id):ep(pat_id));
        
        indivRM = run_mat(sp(pat_id):ep(pat_id),:);
    
        m0 = strfind(indivMenses',[0 1]);
        mf = strfind(indivMenses',[1 0]);
        
        if length(m0) > length(mf)
            mf = [mf length(indivMenses)];
        elseif length(m0) < length(mf)
            m0 = [1 m0];
        end
        
        mdur = mf - m0 + 1;
        valep = mf + tp < length(indivMenses);
        valsp = m0 - tp > 0;
        pullidx = find(mdur > dur_req & valep & valsp);
        
        for j = 1:length(pullidx)
            id = pullidx(j);
            msp = m0(id);
            mep = mf(id);
            preRunMat = indivRM(msp-tp:msp-1,:);
            postRunMat = indivRM(mep+1:mep+tp,:);
            durRunMat = indivRM(msp:mep,:);
            allpreRM(c,:,:) = preRunMat;
            allpostRM(c,:,:) = postRunMat;
            alldurRM(c,:,:) = {durRunMat};
            c = c + 1;
        end
        
        pat_rec(c) = pat_id;
        
    end
            
    % GET MAX CHANGE POINTS DURING (Top 4 max values)
    clear allmaxdur
    for id = 1:size(alldurRM,1)
        tmp = alldurRM{id};
        [v,i] = maxk(tmp(:,1),4);
        newid = sort(i);
        allmaxdur(id,:,:) = tmp(newid,:);
    end

end

%% 
function [switch_ids] = plot_clinical_menses(allpreRM,allpostRM,alldurRM,allmaxdur,tp,th,ci_th)
% Inputs: 
%   * alldurRM: run_mat during menses (nAB, Li and oLB relative abundances)
%   * allpostRM: run_mat post menses (nAB, Li and oLB relative abundances)
%   * allpreRM: run_mat pre menses (nAB, Li and oLB relative abundances)
%   * allmaxdur: maximum values of nAB during menses (nAB, Li and oLB relative abundances)
%   * tp: time period to assess before and after menses
%   * th: threshold used to determine dominance (ex: 0.5 or 0.6)
%   * ci_th: confidence interval thresohld (ex: 0.95 for 95% CI)
%
% Outputs:
%   * switch_ids: indexes of individuals that switched to nAB dominance
%       after menses

    sp_cols = [0.9290 0.6940 0.1250;
            0.5 0.5 0.5;
            0.3010 0.7450 0.9330];

    m0 = tp+1;
    mf = tp+4;
    totL = tp+4+tp;
    
    % menses #1
    switch_ids = NaN(length(alldurRM),1); %-1 = nAB dom before, % 0 = no swittch, % 1 = temp, % 2 = sustained
    trajs = NaN(length(alldurRM),totL,3);
    avgs = NaN(length(alldurRM),3,3);
    for men_id = 1:length(alldurRM)
        tmp_traj = [squeeze(allpreRM(men_id,:,:));
            squeeze(allmaxdur(men_id,:,:));
            squeeze(allpostRM(men_id,:,:))];
        
        tmp_avg = [nanmean(tmp_traj(1:m0-1,:));
            nanmean(tmp_traj(m0:mf,:));
            nanmean(tmp_traj(mf+1:end,:))];
        
        if tmp_avg(1,1) > th % not LB dominant before
            switch_id = 0;
        elseif tmp_avg(1,1) < th && tmp_avg(2,1) > th && tmp_avg(3,1) < th % temp
            switch_id = 1;
        elseif tmp_avg(1,1) < th && tmp_avg(3,1) > th % permenant
            switch_id = 1; 
        elseif tmp_avg(1,1) < th && tmp_avg(2,1) < th && tmp_avg(3,1) < th % no change
            switch_id = 2;
        else
            switch_id = NaN;
        end
        
        switch_ids(men_id) = switch_id;
        trajs(men_id,:,:) = tmp_traj;
        avgs(men_id,:,:) = tmp_avg;
    end
    
    
    tabulate(switch_ids)
    
    resp_names = {'Average','Sensitive','Resilient'};
    
    for resp_id = 1:3
        id_avg = switch_ids == 1 | switch_ids == 2;
        if resp_id == 1
            sz = sum(id_avg);
            avg_traj = squeeze(nanmean(trajs(id_avg,:,:),1));
            std_traj = squeeze(nanstd(trajs(id_avg,:,:),[],1));
        else
            sz = sum(switch_ids == resp_id-1);
            avg_traj = squeeze(nanmean(trajs(switch_ids == resp_id-1,:,:),1));
            std_traj = squeeze(nanstd(trajs(switch_ids == resp_id-1,:,:),[],1));
        end
    
        rm_nan = sum(isnan(avg_traj),2) == 0;
        
        comb_CI = tinv(ci_th,sz-1)*std_traj/sqrt(sz);
        upper_CI = avg_traj + comb_CI;
        lower_CI = avg_traj - comb_CI;
        
        subplot(1,3,resp_id)
        terr = 1:totL;
        terr = terr(rm_nan);
        for i = 1:3
            curve1 = upper_CI(rm_nan,i)';
            curve2 = lower_CI(rm_nan,i)';
            inBetweenRegionX = [terr, fliplr(terr)];
            inBetweenRegionY = [curve1, fliplr(curve2)];
            fill(inBetweenRegionX, inBetweenRegionY, sp_cols(i,:),...
                'FaceAlpha',0.2,'EdgeColor',sp_cols(i,:));
            hold on
            p = plot(terr,avg_traj(rm_nan,i),'LineWidth',1.5,'Color',sp_cols(i,:),...
                'Marker','.','MarkerSize',20);
        end
        patch([m0 mf mf m0], [0 0 1 1], [1 0 0], 'FaceAlpha', 0.1)
        xticks([1:totL])
        xticklabels({'-5','-4','-3','-2','-1','m1','m2','m3','m4','+1','+2','+3','+4','+5'})
        xlabel('Menses Day')
        ylim([0 1])
        set(gca,'fontsize',14)
        ylabel('Abundance')
        title(strcat(resp_names(resp_id),": d0 "  ,num2str(round(sz/sum(id_avg)*100,1)),"% (n = ",num2str(sz),")"),...
            'fontsize',6)
        set(gca,'fontname','Arial') 
        ax = gca;
        ax.LineWidth = 1;
        ax.XColor = 'k'; %
        ax.YColor = 'k'; % 
        ax.FontSize = 7;
        
    end
end