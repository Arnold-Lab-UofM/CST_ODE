%% G3_ABX_Compile_Results_TimeSeries
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Goal: Visualize results for G3_ABX_Run_TimeSeries_Analsis.m and compare
% to the CONRAD BV cohort data published in Gustin et al. 2022: 
% https://pubmed.ncbi.nlm.nih.gov/34560047/
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% REQUIRES: 
%   * Folder and workspace names generated in
%       G3_ABX_Run_TimeSeries_Analysis.m
%   * CONRAD BV clinical data (hard-coded)
%
% OUTPUT:
%   * Figures with the ABX simulations stratefied by BV Clearance types as
%       defined in Gustin et al. 2022
%           - No BV Clearance (no switch to oLB/Li dominance)
%           - Temporary BV Clearance (swith to oLB/Li dominance by end of
%                regimen, return to nAB dominance by day 30)
%           - Sustained BV Clearance (oLB/Li dominance by end of regimen
%                and 30 days after)
%           - Delayed BV Clearance (nAB domninance at end of regimen,
%                oLB/Li dominance 30 days after)
%   * Figure comparing clinical data frequencies to model frequencies
%   * Figure with volcano plot comparing successful vs failed treatment as
%       evaluated by day 30 post-treatment

%% 1. Load the Simulation Time Series data

fd_id = 1; % manuscript has only 1 folder to read in
fdr_names = {'ABXTimeSeries/nAB-HMP_-602lhs-day7-29-Mar-2023'};
ws_nm = '29-Mar-2023-1param-mod-0-for-7d-run.mat';
loc_name = strcat(fdr_names{fd_id},'/',ws_nm);

th = 0.5; % threshold for composition to switch
sp_idx = 1; % Index to swtich to (nAB for menses)
dose_ids = [4]; % index of menses magnitude to plot (4: -2.64, calculated from clinical data)
xlimit = 30; % time period to plot post-abx

%% 2. PLOT TIME SERIES SIMULATION
[f1,f2,f3,f4,fv,Dose_Counts] = plot_multi_panel_ABX_v2(loc_name,th,sp_idx,dose_ids,xlimit);

for k = [f1 f2 f3 f4]
    set(k,'units','inches','position',[1 1 2.5 2])
end
set(fv,'units','inches','position',[1 1 2.5 2])
%% 3. PLOT COMPARISON TO MODEL SIMS

Counts = [Dose_Counts; % Model% HMP
    5	11	10	2]; % CONRAD BV

f5 = figure;
[all_comparison] = compare_abx_sims_to_clinical(Counts);
set(f5,'units','inches','position',[1 1 8 2])

%% 4. COMPILE INTO FINAL FIGURE
fig_list = [f1,f2,f3,f4,f5,fv];
idx = {1,2,3,4,[5 6 7],8};
nrows = 2;
ncols = 4;

FF = combine_plots_to_subplot(fig_list,idx,nrows,ncols);


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ABX SPECIFIC FUNCTIONS 
%   * compare_abx_sims_to_clinical
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [all_comparison] = compare_abx_sims_to_clinical(Counts)
    Frequency = Counts ./ sum(Counts,2);
    n = numel(Counts);
    gr = size(Counts,1);
    spcr = size(Counts,2);
    all_comparison = NaN(gr,spcr);
    superBarP = NaN(n,n);
    for resp_id = 1:spcr
        OBS = [Counts(:,resp_id) sum(Counts,2) - Counts(:,resp_id)];
        pvals = NaN(size(Counts,1),size(Counts,1));
        for c_id1 = 1:gr
            for c_id2 = 1:gr
                observed = [OBS(c_id1,:)', OBS(c_id2,:)'];
                [p,~,~] = chigof(observed);
                pvals(c_id1,c_id2) = p;
            end
        end
    
        pidx = triu(ones(size(Counts,1),size(Counts,1)),1);
    
        tmpp = pvals(pidx == 1);
        all_comparison(:,resp_id) = tmpp;
    
    end
    
    Y = Frequency'*100;
    
    C = [138 58 0;
        152 113 43;
        0 25 127;
        0 7 255]./255;
    
    b = superbar(Y, 'P', superBarP,'BarFaceColor',C);
    tick_loc = [0.75 1.25 1.75 2.25 2.75 3.25 3.75 4.25];
    xticks(tick_loc)
    lbs = {'Model','CONRAD'};
    xticklabels(repmat(lbs,1,4))
    xtickangle(45)
    
    X = reshape(Y',[],1);
    for i = 1:length(X)
        text(tick_loc(i),X(i),num2str(X(i),'%1.f%%'),'VerticalAlignment','top', 'HorizontalAlignment', 'center',...
                'fontname','Arial','Color','w','FontSize',7)
    end
    
    ax = gca;
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.LineWidth = 1;
    ax.FontSize = 7;
    ylabel('Percent Samples')
end


end