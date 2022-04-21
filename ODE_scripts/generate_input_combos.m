%% generate_input_combos.m
%
% INPUT an Mx1 Cell of inputs
% Ex: vectorCell = {[0.18 0.15 0.125 0.1], [-0.05 -0.04 -0.02]};
%
% OUTPUT: Matrix of all parameter value combinations
%
% NOTE: Must be in order in which parameters are listed in model.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [newValueMat] = generate_input_combos(vectorCell)
    
    M = length(vectorCell);
    
    if M == 1
        newValueMat = vectorCell{1};
        if size(newValueMat,1) == 1 % check that is is in column format
            newValueMat = newValueMat';
        end
    elseif M == 2
        p1 = vectorCell{1};
        p2 = vectorCell{2};
        newValueMat = [];
        c = 1;
        for i = 1:length(p1)
            for j = 1:length(p2)
                    newValueMat(c,:) = [p1(i) p2(j)];
                    c = c + 1;  
            end
        end
    elseif M == 3
        p1 = vectorCell{1};
        p2 = vectorCell{2};
        p3 = vectorCell{3};
        newValueMat = [];
        c = 1;
        for i = 1:length(p1)
            for j = 1:length(p2)
                for k = 1:length(p3)
                    newValueMat(c,:) = [p1(i) p2(j) p3(k)];
                    c = c + 1;
                end 
            end
        end
    elseif M == 4
        p1 = vectorCell{1};
        p2 = vectorCell{2};
        p3 = vectorCell{3};
        p4 = vectorCell{4};
        newValueMat = [];
        c = 1;
        for i = 1:length(p1)
            for j = 1:length(p2)
                for k = 1:length(p3)
                    for m = 1:length(p4)
                        newValueMat(c,:) = [p1(i) p2(j) p3(k) p4(m)];
                        c = c + 1;
                    end
                end
            end
        end
    else
        disp([num2str(M) ' is not supported, please enter 1-4 variables'])
    end
end