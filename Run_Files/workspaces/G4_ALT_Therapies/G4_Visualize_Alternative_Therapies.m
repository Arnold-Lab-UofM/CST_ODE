%% G4_Visualization_Alternative_Therapies
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Goal: Compile results from Run_Dose_Duration... and Run_ABX_Prebiotic ...
% into one figure.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% 0. Make Figure Output Directory
figure_output_folder = 'Figures/';
mkdir(figure_output_folder)
set(0, 'DefaultAxesFontName', 'Arial');
%% 1. DOSE-DURATION PLOTTING

result_loc = 'DoseDur/'; % location of result folders

fdr_contents = dir(result_loc);
fdr_full_list = {fdr_contents.name};
fdr_list = fdr_full_list(contains(fdr_full_list,'DoseDur'));
duration = [1 7 14 21 60]; % duration (days)
doses = [-2.64 -2 -1 -0.5 -0.5 -0.125 0]; % dose
th = 0.5; % switch threshold

for i = 1:length(fdr_list)
    fldrnm = fdr_list{i};
    for ev_point = [0 30]
        figtit = strcat(extractBefore(fldrnm,'_'),'-EvalP',num2str(ev_point),...
        '-thr',strrep(num2str(th),'.',''),'.fig');
        f = figure;
        [percent_switch,total_runs,...
            hm_handle] = plot_DoseRegimens_Heatmap_v2(strcat(result_loc,fldrnm),...
            doses,duration,ev_point,th);
        savefig(f,strcat(figure_output_folder,figtit))
        close(f)
    end
end

%% 2. PRE-ABX PLOTTING

result_loc = 'PreAbx/'; % location of result folders

fdr_contents = dir(result_loc);
fdr_full_list = {fdr_contents.name};
fdr_list = fdr_full_list(contains(fdr_full_list,'PreAbx'));
th = 0.5;

for i = 1:length(fdr_list)
    fldrnm = fdr_list{i};
    for ev_point = [0 30]
        figtit = strcat(extractBefore(fldrnm,'_'),'-EvalP',num2str(ev_point),...
        '-thr',strrep(num2str(th),'.',''),'.fig');
        f = figure;
        [results,total_runs,param_values1,...
            param_values2] = plot_Combination_Heatmap(strcat(result_loc,fldrnm),...
            ev_point,th);
        savefig(f,strcat(figure_output_folder,figtit))
        close(f)
    end
end

%% 3. COMPILE INTO MANUSCRIPT FIGURES

eval_point = 'P30'; % select evaluation point (other options is P0)

fig_folder = 'Figures/';
listing = dir(fig_folder);
file_nms = {listing.name};
fig_nms = file_nms(contains(file_nms,'.fig'));

day30 = contains(fig_nms,eval_point); % pulls evaluation point
simType = contains(fig_nms,'Pre'); % keeps track of what type off simulation
fig_order = {'1SS','oLB','Li'}; % order to plot results in

sel_figs_row1 = fig_nms(day30 & simType)';

idxs1 = NaN(length(sel_figs_row1),1); % get resorted indexes for fig_order
for i = 1:length(sel_figs_row1)
    idxs1(i) = find(contains(sel_figs_row1,fig_order(i)));
end

sel_figs_row2 = fig_nms(day30 & ~simType)'; % get resorted indexes for fig_order
idxs2 = NaN(length(sel_figs_row2),1);
for i = 1:length(sel_figs_row2)
    idxs2(i) = find(contains(sel_figs_row2,fig_order(i)))+length(sel_figs_row1);
end

% Compile information to plot into a subplot
fig_list = strcat(fig_folder,[sel_figs_row1; sel_figs_row2]);
idx = num2cell([idxs1;idxs2]);
nrows = 2;
ncols = 3;

finalfigure = combine_plots_to_subplot(fig_list,idx,nrows,ncols);
set(0, 'DefaultAxesFontName', 'Arial');
savefig(finalfigure,strcat('compiled-',eval_point,'-results.fig'))




