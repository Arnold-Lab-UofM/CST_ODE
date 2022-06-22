function xl = sort_run_mat(run_mat)
    mm = NaN(size(run_mat,1),1);
    for i = 1:size(run_mat,1)
        mm(i) = find(run_mat(i,:) == max(run_mat(i,:)));
    end

    [q,idd] = sort(mm,'ascend');

    xl = run_mat(idd,:);
end
 