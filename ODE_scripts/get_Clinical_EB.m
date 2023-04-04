 function [allAssigned,summaryData] = get_Clinical_EB(linearTransMap,SSref,CentroidNames)
   
    SSref = SSref.*8/9;
    
    d = [1,5,9];
    X = linearTransMap(:,d);
    
    for i = 1:size(X,1)
        tmpv = X(i,:);
        distv = sum((SSref - tmpv).^2,2);
        [~,AssignedCentroid] = min(distv);
        allAssigned(i) = AssignedCentroid;
    end
    
    summaryData = array2table(tabulate(allAssigned));
    summaryData.Properties.RowNames = CentroidNames;
    summaryData.Properties.VariableNames = {'Num Identifier','Count','Percentage'};
    disp(summaryData)
    
end