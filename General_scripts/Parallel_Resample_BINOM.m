function [T_sorted_1SE, T_sorted_min, fc_count] = Parallel_Resample_BINOM(num_sets,xblock,yblock,xnames,num_sel,sel_alpha,filename)


% Class 1 and class 2
C1_idx = yblock(:,1) == 1;
C2_idx = yblock(:,2) == 1;

% Selects Sample Size by choosing the size of the smaller class
if sum(C1_idx) == sum(C2_idx)
    ssize = sum(C1_idx);
else % Fix this in a bit
    [min_size, idx] = min([sum(C1_idx) sum(C2_idx)]);
    ssize = min_size;
end




% k-fold is adjusted by sample size
kfold_val = floor(ssize/5); % Number of partitions
inner_kf = floor(ssize/5); % Inner K-fold validation

% Inputs for lassoglm()
NL = 100; % Number of lambdas to try
alpha_val = sel_alpha; % alpha value

seed_val = 1:ceil(num_sets/kfold_val); % How many random data sets to generate


% Split observations for the classes and resample at the same size (size 
% of smaller class)

C1_dat = [xblock(C1_idx,:)];
C2_dat = [xblock(C2_idx,:)];

all_C2_dat = zeros(length(seed_val),ssize,size(xblock,2));
all_C1_dat = zeros(length(seed_val),ssize,size(xblock,2));

for i = 1:length(seed_val)
    rng(seed_val(i));
    [C1,C1_idx] = datasample(C1_dat,ssize,'Replace',false);
    [C2,C2_idx] = datasample(C2_dat,ssize,'Replace',false);
    
    all_C2_dat(i,:,:) = C2;
    all_C1_dat(i,:,:) = C1;
end

% run all combinations of resampled data on lassoglm()

myPool = parpool('local',6);

all_B_1SE = zeros(length(seed_val),kfold_val,size(xblock,2));
all_B_min = zeros(length(seed_val),kfold_val,size(xblock,2));

fc_count = [0 0 0 0]; % Count warning labels 

parfor set_id = 1:length(seed_val)
    ex_C2_dat = squeeze(all_C2_dat(set_id,:,:));
    ex_C1_dat = squeeze(all_C1_dat(set_id,:,:));

    combine_x = [ex_C1_dat;ex_C2_dat];
    combine_y = [ones(ssize,1);zeros(ssize,1)];
    
    rng(set_id) 
    
    c = cvpartition(combine_y,'Kfold',kfold_val);
    for part_id = 1:kfold_val
        temp_idx = training(c,part_id);
        xblock = combine_x(temp_idx,:);
        yblock = combine_y(temp_idx,:);
        [fc,B_1SE,B_min,~,~] = LASSO_ELASTIC_BINOMIAL(xblock,yblock,xnames,inner_kf,NL,alpha_val);
        all_B_1SE(set_id,part_id,:) = B_1SE;
        all_B_min(set_id,part_id,:) = B_min;
        
        fc_count = fc_count + fc;
        
        [set_id part_id];
    end
end

delete(myPool)

% Reshape workspace and evaluate

c = 1;
re_B_1SE = zeros(size(all_B_1SE,1)*size(all_B_1SE,2),size(all_B_1SE,3));
re_B_min = zeros(size(all_B_1SE,1)*size(all_B_1SE,2),size(all_B_1SE,3));
for i = 1:size(all_B_1SE,1)
    for j = 1:size(all_B_1SE,2)
        re_B_1SE(c,:) = all_B_1SE(i,j,:);
        re_B_min(c,:) = all_B_min(i,j,:);
        c = c + 1;
    end
end

fr_sel_1SE = sum(re_B_1SE ~= 0)./(size(all_B_1SE,1)*size(all_B_1SE,2));
fr_sel_min = sum(re_B_min ~= 0)./(size(all_B_1SE,1)*size(all_B_1SE,2));


% Compile Data and Save
% load('mod_names.mat')
% xnames = namfd;
T = array2table([fr_sel_1SE', fr_sel_min'],'RowNames',cellstr(xnames), ...
    'VariableNames',{'EN_1SE','EN_min'});

t1 = T(T.EN_1SE >= 0.7,1);
t2 = T(T.EN_min >= 0.7,1);

T_sorted_1SE = sortrows(T,1,'descend');
T_sorted_min = sortrows(T,2,'descend');

tot_mdls = size(all_B_1SE,1)*size(all_B_1SE,2);

ttb = strcat('Resampled_EN_summary_',num2str(tot_mdls),'Mdls_',extractBefore(filename,'.mat'),'.xlsx');
writetable(T_sorted_1SE,ttb,'WriteRowNames',true)
save(strcat(extractBefore(ttb,'.xlsx'),'.mat'))

%% Plot Result
% num_sel = 15;
[newb idxF] = maxk(T.EN_1SE,num_sel);

new_var = exp(re_B_1SE(:,idxF));
freq_var = T.EN_1SE(idxF);

figure(1)
subplot(1,2,1)
boxplot(new_var, 'orientation', 'horizontal');
set(gca,'ytick',1:length(idxF),'yticklabels',xnames(idxF));
hold on
xline(1,'k')
xlabel('Odds Ratio')
title(strcat("Frequency for top ", num2str(num_sel),  " Features in ", num2str(tot_mdls)," 1SE Mdls"))
set(gca,'xscale','log')

subplot(1,2,2)
b = barh(freq_var);
b.FaceColor = [0.75 0.75 0.75];
set(gca,'ytick',1:length(idxF),'yticklabels',xnames(idxF));
hold on
xline(0.7,'k')
xlabel('Frequency of Selection')
title(strcat("Frequency for top ", num2str(num_sel),  " Features in ", num2str(tot_mdls)," 1SE Mdls"))
savefig(gcf,strcat(extractBefore(filename,'.mat'),'-1SE.fig'))

close
[newb idxF] = maxk(T.EN_min,num_sel);

new_var = exp(re_B_min(:,idxF));
freq_var = T.EN_min(idxF);

figure(2)
subplot(1,2,1)
boxplot(new_var, 'orientation', 'horizontal');
set(gca,'ytick',1:length(idxF),'yticklabels',xnames(idxF));
hold on
xline(1,'k')
xlabel('Odds Ratio')
title(strcat("Frequency for top ", num2str(num_sel),  " Features in ", num2str(tot_mdls)," Min MSE Mdls"))
set(gca,'xscale','log')

subplot(1,2,2)
b = barh(freq_var);
b.FaceColor = [0.75 0.75 0.75];
set(gca,'ytick',1:length(idxF),'yticklabels',xnames(idxF));
hold on
xline(0.7,'k')
xlabel('Frequency of Selection')
title(strcat("Frequency for top ", num2str(num_sel),  " Features in ", num2str(tot_mdls)," Min MSE Mdls"))
savefig(gcf,strcat(extractBefore(filename,'.mat'),'-MINSE.fig'))

%https://www.mathworks.com/help/stats/datasample.html
end