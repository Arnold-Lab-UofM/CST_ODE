%% get_SS_info_3sp.m
% [poss_SS, mat, poss_SSnames, mat_names,mat_num] = get_SS_info(SS_vector,plotFlag)
%
% GOAL: Assess the stability of each input parameter set and list if
%   multiple states are stable ("steady state configurations").
%
% INPUTS:
%   * SS_vector: cell vector of stable steady states for multiple parameter 
%       sets
%   * plotFlag: plot summary (plots if true)
%
% OUTPUTS:
%   * poss_SS: matrix of 0's and 1's for what steady states are stable for
%       a given parameter set [N x 8], N = number of parameter sets
%   * mat: matrix of 0's and 1's for what configuration is assigned to a
%       certain parameter set [N x 20], 20 = number of unique
%       configurations of the 8 steady states
%   * poss_SSnames: names of the 8 steady states
%   * mat_names: names of the 20 configurations
%   * mat_num: numerical value assigned to each of the 20 configurations
%   * mat_code: association with healthy or BV state
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% Jan 12, 2020
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [poss_SS, mat, poss_SSnames, mat_names,mat_num,mat_code] = get_SS_info_3sp(SS_vector,plotFlag)
%%
    poss_SSnames = {'Co-elimination','Only BV','Only Li','Only Lc',...
        'Bv and Li', 'BV and Lc', 'Li and Lc', 'Co-existence'};
    poss_SS = zeros(size(SS_vector,1),8);
    for i = 1:size(SS_vector,1)
        temp = SS_vector{i}{1};
        temp1s = temp > 0;
        for j = 1:size(temp1s)
            tp = temp1s(j,:);
            if sum(tp) == 0
                poss_SS(i,1) = 1;
            elseif sum(tp) == 1
                idx = find(tp,1);
                poss_SS(i,idx+1) = 1;
            elseif sum(tp) == 2
                if sum(tp(1:2)) == 2
                    poss_SS(i,5) = 1;
                elseif sum(tp([1,3])) == 2
                    poss_SS(i,6) = 1;
                else
                    poss_SS(i,7) = 1;
                end
            else
                poss_SS(i,8) = 1;
            end
        end
    end

    %% %% Common Patterns
    idx234 = poss_SS(:,2) & poss_SS(:,3) & poss_SS(:,4);
    idx27 = poss_SS(:,2) & poss_SS(:,7);
    idx57 = poss_SS(:,5) & poss_SS(:,7);
    idx45 = poss_SS(:,4) & poss_SS(:,5);
    idx36 = poss_SS(:,3) & poss_SS(:,6);
    idx56 = poss_SS(:,5) & poss_SS(:,6);
    idx67 = poss_SS(:,6) & poss_SS(:,7);
    idx7o = poss_SS(:,7) & sum(poss_SS,2) == 1;
    idx6o = poss_SS(:,6) & sum(poss_SS,2) == 1;
    idx5o = poss_SS(:,5) & sum(poss_SS,2) == 1;
    idx4o = poss_SS(:,4) & sum(poss_SS,2) == 1;
    idx3o = poss_SS(:,3) & sum(poss_SS,2) == 1;
    idx2o = poss_SS(:,2) & sum(poss_SS,2) == 1;
    idx8o = poss_SS(:,8) & sum(poss_SS,2) == 1;
    idx23 = poss_SS(:,2) & poss_SS(:,3) & sum(poss_SS,2) == 2;
    idx34 = poss_SS(:,3) & poss_SS(:,4) & sum(poss_SS,2) == 2;
    idx24 = poss_SS(:,2) & poss_SS(:,4) & sum(poss_SS,2) == 2;
    idx28 = poss_SS(:,2) & poss_SS(:,8) & sum(poss_SS,2) == 2;
    idx38 = poss_SS(:,3) & poss_SS(:,8) & sum(poss_SS,2) == 2;
    idx48 = poss_SS(:,8) & poss_SS(:,4) & sum(poss_SS,2) == 2;

%     mat = [idx234,idx27,idx36,idx45,idx56,idx57,idx67,idx2o,idx3o,idx4o,idx5o,idx6o,idx7o,idx8o];
    mat = [idx234,idx23,idx24,idx34,idx2o,idx3o,idx4o,idx27,idx36,idx45,...
        idx57,idx56,idx67,idx5o,idx6o,idx7o,idx28,idx38,idx48,idx8o];
%     sum(mat,'all')

    mat_names = {'3SS - solos', '2SS - Bv or Li', '2SS - Bv or Lc', ...
        '2SS - Li or Lc', '1SS - Bv', '1SS - Li', '1SS - Lc',...
        '2SS - Bv or LiLc', '2SS - Li or BvLc', '2SS - Lc or BvLi', ...
        '2SS - BvLi or LiLc', '2SS - BvLi or BvLc', '2SS - BvLc or LiLc',...
        '1SS - BvLi', '1SS - BvLc', '1SS - LiLc', '2SS - Bv or BvLiLc',...
        '2SS - Li or BvLiLc', '2SS - Lc or BvLiLc','1SS - BvLiLc'};
    
    mat_code = {'ANY', 'ANY','ANY','HEALTHY','BV','HEALTHY','HEALTHY',...
        'ANY','ANY','ANY','ANY','ANY','ANY','BV','BV','HEALTHY','ANY',...
        'ANY','ANY','ANY'};
    
    for i = 1:size(mat,1)
        tmp = find(mat(i,:));
        if isempty(tmp)
            mat_num(i) = NaN;
        elseif length(tmp) > 1
            mat_num(i) = NaN;
            disp('Run does not fit one category')
            disp(i)
        else
            mat_num(i) = tmp;
        end
    end
    
    if plotFlag
        figure
        bar(sum(mat)./sum(mat,'all'))
        hold on
        text(1:size(mat,2),sum(mat)./sum(mat,'all'),string(round(sum(mat)./sum(mat,'all'),3)'),'vert','bottom','horiz','center'); 
        xticks(1:size(mat,2))
        xticklabels(mat_names)
        xtickangle(90)
        ylabel('Frequency')
        title('Combinations of Possible Steady-States')
        set(gca,'fontsize',12)
    end
end
