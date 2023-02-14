%% [PrismFormat,SummaryStats] = plot_Volcano(matrix1,matrix2,alpha,offset,param_names,classes)
%
% INPUT:
%   * matrix1 & matrix2: parameter sets you want to compare
%   * alpha: significance threshold
%   * offset: plotting of data labels distance from point
%   * classes: names of the groups being compared
%
% OUTPUT:
%   * PrismFormat: matrix to copy and paste into prism to easily plot with
%       color-coded points (by statistical significance, positive or negative)
%   % SummaryStates: Summary statistics between the two groups
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% Jan 20, 2023 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [PrismFormat,SummaryStats] = plot_Volcano(matrix1,matrix2,alpha,offset,param_names,classes)
    
    Ns =[size(matrix1,1) size(matrix2,1)];

    % Compare each parameter between the two groups
    for k = 1:length(param_names)
        tmp1 = matrix1(:,k);
        tmp2 = matrix2(:,k);
        
        p = ranksum(tmp1,tmp2); % Non-parametric comparison
        svp(k) = p;
        tmp = NaN(max(Ns),2);
        tmp(1:Ns(1),1) = tmp1;
        tmp(1:Ns(2),2) = tmp2;

        X(:,k) = nanmean(tmp)';
        ylabel(param_names(k))
    
        [mean1,mean2]=  mean_rank(tmp1,tmp2);
    
        DIF(k) = mean2 - mean1; % compare mean rank
    end
    
    SummaryStats = rows2vars(array2table(X,'RowNames',classes,'VariableNames',...
        param_names));
    
    % ~~~~~~~~~~ VOLCANO ~~~~~~~~~~
    [FDR] = mafdr(svp,'BHFDR',true); % FDR multiple comparison p-value adjustment
    
    idxLR = DIF < 0;
    idxsig = FDR < alpha;
    
    idxR = ~idxLR & idxsig;
    idxL = idxLR & idxsig;
    
    plot(DIF(idxL),-log10(FDR(idxL)),'o','MarkerFaceColor',[59 110 178]/255,'MarkerEdgeColor','k',...
        'MarkerSize',10)
    hold on
    plot(DIF(idxR),-log10(FDR(idxR)),'o','MarkerFaceColor',[194 171 131]/255,'MarkerEdgeColor','k',...
        'MarkerSize',10)
    hold on
    plot(DIF(~idxsig),-log10(FDR(~idxsig)),'o','MarkerFaceColor',[200 200 200]/255,'MarkerEdgeColor','k',...
        'MarkerSize',10)
    yline(-log10(alpha),'LineStyle',':')
    xline(0,'LineStyle',':','LineStyle',':')
    text(DIF+offset,-log10(FDR)+offset,param_names)
    xlabel('Rank Difference')
    ylabel('-log10(q-value)')
    set(gca,'fontsize',14,'fontname','arial')
    title(classes)
    
    %
    PrismFormat = NaN(length(param_names),4);
    
    PrismFormat(:,1) = DIF;
    PrismFormat(idxL,2) = -log10(FDR(idxL));
    PrismFormat(idxR,4) = -log10(FDR(idxR));
    PrismFormat(~idxsig,3) = -log10(FDR(~idxsig));
end