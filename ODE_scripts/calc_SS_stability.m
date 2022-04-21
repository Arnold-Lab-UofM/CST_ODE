%% calc_SS_stability.m
% Use this to calculate stability of 'GUASE'-type or 'CAlpha' type
% generalized Lotka Volterra Equations
% [StableStates,SSval,eigval] = calc_SS_stability(num_sp,params,S,Jmat,Type)
% Input:
% num_sp = numbers of species
% params = parameter set (growth rates, carrying capacity, interaction
%   terms)
% S = steady-state solutions calculated by symbolic_solns.m
% Jmat = Jacobian from symbolic solns.m
% Type = 'GAUSE' or 'CAlpha' dependence on type of gLV formulation
%
% Output: 
% StableStates = all biologically relevant steady states
% SSval = all steady states (unstable and negative)
% eival = eigen values for each steady-state
%
% Questions can be directed at Christina Lee (chyylee@umich.edu)
% 2/19/2021
function [StableStates,SSval,eigval,UnstableStates] = calc_SS_stability(num_sp,params,S,Jmat,Type)
     N = num_sp;

    %% Example 3 species
    sstates = matlabFunction(struct2array(S));
    Jm = matlabFunction(Jmat);

    %% Reformat parameter order
    if contains(lower(Type),'gause')
      % Reformat parameter order
        p_grow = params(1:N);
        K = params(N+1:2*N);
        fx = reshape(params(2*N+1:end),[N N]);
        p_int = zeros(size(fx));
      if contains(lower(Type),'fx')
            for i = 1:N
                for j = 1:N
                    p_int(i,j) = (K(i) - fx(i,j)*K(i))/(fx(j,i)*K(j));
                end
            end
      else
          p_int = fx;
      end
        
        % Remove dummy self-interaction terms
        idx_rm = linspace(1,N^2,N);
        newbeta = reshape(p_int,[],1);
        newbeta(idx_rm) = [];

        ssparam = num2cell(horzcat(newbeta',K,p_grow));

        SSval = arrayfun(sstates,ssparam{1:end-N},'UniformOutput',false);
    elseif contains(lower(Type),'calpha')
        p_grow = params(1:N);
        K = params(N+1:2*N);
        fx = reshape(params(2*N+1:end),[N N]);
        p_int = zeros(size(fx));
        for i = 1:N
            for j = 1:N
                if i~= j
                    p_int(i,j) = fx(i,j);
                else
                    p_int(i,j) = - p_grow(i)/K(i);
                end
            end
        end
            idx_rm = linspace(1,N^2,N);
            newbeta = reshape(p_int,[],1);
            newbeta(idx_rm) = [];

            ssparam = num2cell(horzcat(newbeta',K,p_grow));

            SSval = arrayfun(sstates,ssparam{:},'UniformOutput',false);
    elseif contains(lower(Type),'ralpha')
            p_grow = params(1:N);
            p_int = params(N+1:end);

%             idx_rm = linspace(1,N^2,N);
            newbeta = reshape(p_int,[],1);
%             newbeta(idx_rm) = [];

            ssparam = num2cell(horzcat(newbeta',p_grow));

            SSval = arrayfun(sstates,ssparam{:},'UniformOutput',false);
        
    else
        disp('Please enter: FXGAUSE, GAUSE, RALPHA or CALPHA for input Type')
    end

    %% Calculate analytical SS values and stability
  
    SSval = SSval{1};
    eigval = NaN(size(SSval));
    % Loop through all Steady-states
    for i = 1:size(SSval,1)
        if sum(isnan(SSval(i,:)))>0
            eigval(i,:) = NaN(size(SSval(i,:)));
        elseif sum(SSval(i,:) == Inf)>0
            eigval(i,:) = NaN(size(SSval(i,:)));
        else
            iv = horzcat(ssparam,num2cell(SSval(i,:)));
            Jval = arrayfun(Jm,iv{:},'UniformOutput',false);
            eigval(i,:) = eig(Jval{1});
        end
    end

    idx = max(real(eigval),[],2)<=0;
    idx_neg = min(SSval,[],2)<0;

    StableStates = SSval(idx & ~idx_neg,:);
    UnstableStates = SSval(~idx & ~idx_neg,:);

end

