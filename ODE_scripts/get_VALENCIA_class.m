function [nmf] = get_VALENCIA_class(run_mat)
    n1 = size(run_mat,1);
    
    if n1 > 1
        run_mat = sort_run_mat(run_mat);
    end
    
    SSnms = {'[NO] CST-IV','[Li] CST-III','[oLB] CST-I/II/V','NoStbSS'};
    
    SSref = [0.911801263	0.059240659	0.028958078;
        0.146456352	0.758938299	0.09460535;
        0.153186234	0.09522381	0.751589956];


    SS_type = ones(size(run_mat,1),1);
    SS_typenms = cell(size(run_mat,1),1);
    SS_map = zeros(size(run_mat,1),size(SSref,1));
    nm = [];
    for i = 1:n1
        tmp = run_mat(i,:);
        if sum(isnan(tmp)) > 0
            id = 4;
        else

            dif = sum((SSref - repmat(tmp,size(SSref,1),1)).^2,2);

            [val,id] = min(dif);
        end

        SS_type(i) = id;
        SS_typenms(i) = SSnms(id);
        SS_map(i,id) = 1;
        
        if i > 1
            nm = strcat(nm," or ", SSnms{id});
        else
            nm = SSnms{id};
        end
    end
    
    nmf = strcat(num2str(n1),"SS: ", nm);

end