
function s = lhs_ode_norm(xmean, xsd, nsample)
%============================================================================
%
% LHS from normal distribution, no correlation
% method of Stein - s = latin_hs(xmean,xsd,nsample,nvar)
%
% Stein, M. 1987. Large Sample Properties of Simulations Using Latin
%                 Hypercube Sampling. Technometrics 29:143-151
%
% Inputs:
%
%	xmean
%
%	xsd
%
%   nsample: sample size
%
% Output:
%
%	s: a vector of size nsample with random values selected, using the LHS
%		 algorithm, from a normal probability distribution.
%
%============================================================================

if (nsample == 1)
    s = xmean;
    return
end

ran = rand(nsample, 1);
s = zeros(nsample, 1);

% method of Stein
idx = randperm(nsample);

% probability of the cdf
prob = (idx'-ran)/nsample;

% this can be replaced by any inverse distribution function
s = xmean + ltqnorm(prob).* xsd;

