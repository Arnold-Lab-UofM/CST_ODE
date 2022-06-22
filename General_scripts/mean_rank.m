    function [mean1,mean2]=  mean_rank(tmp1,tmp2)
        x = [tmp1;tmp2];

        for i = 1:length(tmp1)
            value_to_find = tmp1(i);
            position=find(value_to_find==sort(x),1,'first');
            all_pos1(i) = position;
        end

        for i = 1:length(tmp2)
            value_to_find = tmp2(i);
            position=find(value_to_find==sort(x),1,'first');
            all_pos2(i) = position;
        end

        mean1 = mean(all_pos1);
        mean2 = mean(all_pos2);
    end