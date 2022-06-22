%% LASSO_ELASTIC_BINOMIAL.m


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Lee, Arnold Lab, University of Michigan, Biomedical Engineering
% April 10th, 2020
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [flag_count,B_1SE,B_min,xnames_min,xnames_1SE] = LASSO_ELASTIC_BINOMIAL(xblock,yblock,xnames,kf_val,NL,ALPHA)


%% Call Original LASSO Function
%evrimovepath('bottom') % Moves PLS_Toolbox to bottom of path

% Runs Elastic Net by Setting alpha 0.5, alpha = 1 would be Lasso, alpha =
% 0 would be equivalent to ridge regresssion

% Currently runs 10-fold cross-validation and calculates MSE (CV = 10)

% NumLamda tells the function how many different values of Lambda to try
% (Lambda is a regularization parameter, as it increases the number of
% non-zero components of predictors, beta, decreases)
% kf_val = 10;
% NL = 100;
% ALPHA = 0.8;

% Get Boolean Input 
ybool = yblock(:,1) == 1;

% Mean center and variance scale
Xnorm = zscore(xblock);

[B, FitInfo] = lassoglm(Xnorm,ybool,'binomial','CV',kf_val,'NumLambda',NL,'alpha', ALPHA);
% lassoPlot(B,FitInfo,'PlotType','CV');
% legend('show') % Show legend

% B are the coefficients for the different values of Lambda
% Fitinfo is a structure that stores model information

%% Select Lambda value and associated predictor coefficients
% flag counter:

if FitInfo.IndexMinDeviance == length(FitInfo.Lambda) && FitInfo.Index1SE ~= length(FitInfo.Lambda)
    % If there is an issue, the code manually pulls the index of lambda
    % value that is 1 std deviation from the minimum error (both
    % deviance or standard error metric) and displays error message
    
%     % Based on deviance:
%     deviance_vals = FitInfo.Deviance;
%     std_dev_dev = std(deviance_vals);
%     [minVal, IDX_dev] = min(abs(deviance_vals - (min(deviance_vals) + std_dev_dev)));

    % Based on 1SE:
    se_vals = FitInfo.SE;
    std_dev_se = std(se_vals);
    [minVal, IDX_se] = min(abs(se_vals - (min(se_vals) + std_dev_se)));
    
    B_1SE = B(:,FitInfo.Index1SE);
    B_min = B(:,IDX_se);
    
    % SE associated with selected model
    cvMSE_1SD = FitInfo.SE(IDX_dev);
    cvMSE_min = FitInfo.SE(IDX_se);
    
    % Displays warning flag
    disp('WARNING: the mininmum error model had all zero coefficients')
    flag_count = [0 0 1 0];
elseif FitInfo.Index1SE == length(FitInfo.Lambda) && FitInfo.IndexMinDeviance == length(FitInfo.Lambda)
   % If there is an issue, the code manually pulls the index of lambda
    % value that is 1 std deviation from the minimum error (both
    % deviance or standard error metric) and displays error message
    
    % Based on deviance:
    deviance_vals = FitInfo.Deviance;
    std_dev_dev = std(deviance_vals);
    [minVal, IDX_dev] = min(abs(deviance_vals - (min(deviance_vals) + std_dev_dev)));

    % Based on 1SE:
    se_vals = FitInfo.SE;
    std_dev_se = std(se_vals);
    [minVal, IDX_se] = min(abs(se_vals - (min(se_vals) + std_dev_se)));
    
    B_1SE = B(:,IDX_dev);
    B_min = B(:,IDX_se);
    
    % SE associated with selected model
    cvMSE_1SD = FitInfo.SE(IDX_dev);
    cvMSE_min = FitInfo.SE(IDX_se);
    
    % Displays warning flag
    disp('WARNING: the mininmum error an 1se model had all zero coefficients')
    flag_count = [0 0 0 1];
elseif FitInfo.Index1SE == length(FitInfo.Lambda) && FitInfo.IndexMinDeviance ~= length(FitInfo.Lambda)
    % If lassoglm() has no issues then proceeds as normal
    
    se_vals = FitInfo.SE;
    std_dev_se = std(se_vals);
    [minVal, IDX_se] = min(abs(se_vals - (min(se_vals) + std_dev_se)));
    
    B_1SE = B(:,IDX_se);
    B_min = B(:,FitInfo.IndexMinDeviance);
    
    % SE associated with selected model
    cvMSE_1SD = FitInfo.SE(IDX_se);
    cvMSE_min = FitInfo.SE(FitInfo.IndexMinDeviance);
    disp('WARNING: the 1se model had all zero coefficients')
    flag_count = [0 1 0 0];
else
    B_1SE = B(:,FitInfo.Index1SE);
    B_min = B(:,FitInfo.IndexMinDeviance);
    
    % SE associated with selected model
    cvMSE_1SD = FitInfo.SE(FitInfo.Index1SE);
    cvMSE_min = FitInfo.SE(FitInfo.IndexMinDeviance);
    flag_count = [1 0 0 0];    
end

% B_1SE = B(:,FitInfo.Index1SE);
% B_min = B(:,FitInfo.IndexMinDeviance);

% Cross-validation error
% cvMSE_1SD = FitInfo.SE(FitInfo.Index1SE);
% cvMSE_min = FitInfo.SE(FitInfo.IndexMinDeviance);

% load('mod_names.mat')
% xnames = namfd;
% Selected Parameter Names
xnames_1SE = xnames(B_1SE ~= 0);
xnames_min = xnames(B_min ~= 0);

% Reformatted Matrix with only Selected Predictors
xblock_1SE = xblock(:,B_1SE ~= 0);
xblock_min = xblock(:,B_min ~= 0);



%% Call Model to Make Predictions
% cnst = FitInfo.Intercept(FitInfo.Index1SE);
% B1 = [cnst;B_1SE];
% 
% preds = glmval(B1,xblock,'logit');
% histogram(ybool - preds)
% title('Residuals from lassoglm model')

% If in the future you want to feed other testing data, use these functions
% to make new predictions.



end
%% Visualize Predictors (Bar Graph)

% bar([B_1SE(B_1SE~=0)])
% set(gca,'xtick',1:length(xnames_1SE),'xticklabels',xnames_1SE)
% xtickangle(270)
% xlabel('Predictor')
% ylabel('Coefficient')
% title(['MSE 1SE = ', num2str(round(cvMSE_1SD,3))])
% 
% 
% array2table(B_1SE(B_1SE~=0),'RowNames',cellstr(xnames_1SE),'VariableNames', ...
%     {'Nonzero_x1SE'})


% %% Save workspace
% if ALPHA == 0.5
%     modelfilettl = strcat(datestr(today()),strtok(extractAfter(filettl,11),'.'),'_ELASTICNET_LamMin.mat');
%     xblock = xblock_min;
%     xnames = xnames_min;
%     save(modelfilettl)
%     
%     modelfilettl = strcat(datestr(today()),strtok(extractAfter(filettl,11),'.'),'_ELASTICNET_Lam1SE.mat');
%     xblock = xblock_1SE;
%     xnames = xnames_1SE;
%     save(modelfilettl)
% elseif ALPHA == 1
%     modelfilettl = strcat(datestr(today()),strtok(extractAfter(filettl,11),'.'),'_LASSO_LamMin.mat');
%     xblock = xblock_min;
%     xnames = xnames_min;
%     save(modelfilettl)
%     
%     
%     modelfilettl = strcat(datestr(today()),strtok(extractAfter(filettl,11),'.'),'_LASSO_Lam1SE.mat');
%     xblock = xblock_1SE;
%     xnames = xnames_1SE;
%     save(modelfilettl)
% else
%     modelfilettl = strcat(datestr(today()),strtok(extractAfter(filettl,11),'.'),'_alpha=',num2str(ALPHA),'LamMin.mat');
%     xblock = xblock_min;
%     xnames = xnames_min;
%     save(modelfilettl)
%     
%     modelfilettl = strcat(datestr(today()),strtok(extractAfter(filettl,11),'.'),'_alpha=',num2str(ALPHA),'Lam1SE.mat');
%     xblock = xblock_1SE;
%     xnames = xnames_1SE;
%     save(modelfilettl)
% end

%% Other info/helpful Links
% https://www.mathworks.com/help/stats/lasso-and-elastic-net.html
% https://www.mathworks.com/help/stats/classificationpartitionedlinear.kfoldloss.html
% https://www.mathworks.com/help/stats/crossval.html
% http://wiki.eigenvector.com/index.php?title=Faq_good_way_to_add_or_remove_PLS_Toolbox_from_my_path

