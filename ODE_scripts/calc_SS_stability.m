%% calc_SS_stability.m
% [StableStates,SSval,eigval] = calc_SS_stability(N,params,S,Jmat,Type)
%
% Use this to calculate stability of 'GUASE'-type, 'CAlpha', 'Ralpha' type
% generalized Lotka Volterra Equations
%
% Input:
%   * N = numbers of species
%   * params = parameter set (growth rates, carrying capacity, interaction
%       terms)
%   * S = steady-state solutions calculated by symbolic_solns.m
%   * Jmat = Jacobian from symbolic solns.m
% Output: 
%   * StableStates = all biologically relevant steady states
%   * SSval = all steady states (unstable and negative)
%   * eival = eigen values for each steady-state
%   * UnstableStates = non-biologically relevant steady states
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% Feb 19, 2021
% Update: Jan 20, 2023 (Removed all but gLV option, changed cell handling)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [StableStates,SSval,eigval,UnstableStates] = calc_SS_stability(N,params,S,Jmat)

    %% Example 3 species
    sstates = matlabFunction(struct2array(S));
    Jm = matlabFunction(Jmat);

    %% Reformat parameter order
    p_grow = params(1:N);
    p_int = params(N+1:end);

    newbeta = reshape(p_int,[],1);

    ssparam = num2cell(horzcat(newbeta',p_grow));

    SSval = cell2mat(arrayfun(sstates,ssparam{:},'UniformOutput',false));

    %% Calculate analytical SS values and stability
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
            eigval(i,:) = eig(cell2mat(Jval));
        end
    end

    idx = max(real(eigval),[],2)<=0; % Eigen values less than or equal to zero
    idx_neg = min(SSval,[],2)<0; % remove negative SS options

    StableStates = SSval(idx & ~idx_neg,:);
    UnstableStates = SSval(~idx & ~idx_neg,:);

end


