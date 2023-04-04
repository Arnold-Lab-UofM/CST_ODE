
function rm_indx = remove_outlier_states(sel_nets,ybase,S,Jmat,tspan)
   
%%
    [~,mxidx] = min(ybase);
    tall = cell(size(sel_nets,1),1);
    yall = cell(size(sel_nets,1),1);
    yout = zeros(size(sel_nets,1),3);
    ymax = zeros(size(sel_nets,1),1);
    warnall = strings(size(sel_nets,1),1);
    
    icth = 0.5;
    for net_id = 1:size(sel_nets,1)
        base_params = sel_nets(net_id,:);
        
        % Determine base SS
        [OGSS,~,~] = calc_SS_stability(length(ybase),base_params,S,Jmat);
        [~,midxss] = max(OGSS,[],2);
        [~,id] = intersect(midxss,mxidx);
    
        if isempty(id)
            id = 1;
        end
        
        y0 = OGSS(id,:) + ybase*sum(OGSS(id,:)); % Set initial conditions based on input
    
        [t, y] = ode45(@lhs_ode_gLV,tspan,y0,[],base_params);
    
        yout(net_id,:) = y(end,:)./sum(y(end,:),2);
    
        ybv = max(y(:,1)./sum(y,2));
        ymax(net_id) = ybv;
        [warnmsg1, ~] = lastwarn; 
        warnall(net_id) = warnmsg1;
        lastwarn('')
        
    end
    
    [~,spidx] = min(ybase); % checks for dom species idx (will be negative, thus the min)
    if spidx ~= 1
        idx_rm1 = yout(:,1)>icth; % remove parameter sets that switch without perturbation;
        idx_rm2 = ymax > 0.4;
    else
        idx_rm1 = yout(:,1)<icth; % remove parameter sets that switch without perturbation;
        idx_rm2 = ymax < 0.4;
    end
    idx_rm3 = warnall ~= "";
    
    rm_indx = idx_rm2 | idx_rm1 | idx_rm3;
end