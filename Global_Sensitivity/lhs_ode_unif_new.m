
function s = lhs_ode_unif_newer(xmin, xmax, nsample, logUniform)
%============================================================================
%
% LHS from uniform distribution, no correlation method of Stein -
% s=latin_hs(xmean,xsd,nsample,nvar)
%
% Stein, M. 1987. Large Sample Properties of Simulations Using Latin Hypercube
% Sampling.  Technometrics 29:143-151
%
% NOTE:
%    This does not work properly for log uniform for a range where xmin < 0,
%    xmax > 0 and log(abs(xmin)) < log(xmax).
%
% Inputs:
%
%	xmin: The minimum value for the range.
%
%	xmax: The maximum value for the range.
%
%   nsample: sample size
%
%   logUniform: true to use log uniform, false to use uniform.
%
% Output:
%
%	s: a vector of size nsample with random values selected, using the LHS
%	algorithm, from a uniform or log uniform probability distribution.
%
%============================================================================

if xmin >= xmax
    msgIdent = 'ODE:UNIF:InvalidRange';
    msgfmt = 'Invalid range for getting random uniform or log uniform number,';
    msgfmt = [msgfmt ' xmin, %d, is >= xmax, %d'];
    msg = sprintf(msgfmt, xmin, xmax);
    throw(MException(msgIdent, msg));  
end

if logUniform
    if xmin == 0
        msgIdent = 'ODE:UNIF:InvalidXmin';
        msg = 'Xmin is 0.0 for getting random log uniform number';
        throw(MException(msgIdent, msg));
    end
    
    if xmax == 0
        msgIdent = 'ODE:UNIF:InvalidXmax';
        msg = 'Xmax is 0.0 for getting random log uniform number';
        throw(MException(msgIdent, msg));
    end
end


if (nsample == 1)
    s = xmin;
    return;
end

% Initialize the result vector.
s = zeros(nsample, 1);

% Get a column vector of random probability values (values in [0, 1]).
ran = rand(nsample, 1);

% Get a column vector of sample ordinals in random order.
idx = randperm(nsample)';

% Define a probability for each subinterval of the range.
% For the ith subinterval (the element of idx = i), when the corresponding
% element of ran is 0, prob is (i - 0)/nsample = 1/nsample. When ran is 1
% prob is (i - 1)/nsample. So the range for prob for subinterval i is 
% (i - 1)/nample to i/nsample. Suppose nsample = 2, and idx = [2, 1].
% prob(1) is in the range ((2 - 1)/2, 2/2) = (0.5 to 1) and prob(2) is in
% the range ((1 - 1)/2, 1/2) = (0, 0.5).
%
% This effectively divides the interval (0, 1) into nsample subintervals
% and chooses a value in that subinterval. Randomizing the order of idx
% randomizes the order of those subintervals.
prob = (idx - ran(:))/nsample;

if (~logUniform)
    % Linear scale.
    s(:) = xmin + prob.* (xmax - xmin);
else
    % Log scale

    % If part or all of the range is negative, adjust it so its positive.
    posMin = xmin;
    posMax = xmax;
    bothNegative = false;
    xminNegativeOnly = false;
    if xmin < 0 && xmax < 0
        % The entire range is negative. 
        % Make it all positive by switching the end points, as positive
        % values. This also preserves the relative order of magnitude.
        posMin = abs(xmax);
        posMax = abs(xmin);
        bothNegative = true;
    elseif xmin < 0 && xmax > 0
        % Part of the range is negative and part positive.
        % Make it all positive by flipping it about xmax.
        % xmax is the new min (since it is positive) and the range between the
        % new min and new max is the same as before:
        % posMax - posMin = (xmax + (xmax - xmin)) - xmax = xmax - xmin.
        % This also preserves the relative order of magnitude, except in
        % the case where -1 < xmin < 0 and xmax > 1. Ex. xmin = -0.1 and
        % xmax = 100.
        posMax = xmax + (xmax - xmin);
        posMin = xmax;
        xminNegativeOnly = true;
    end
    
    % Exponents (log of min and max) in the same range: both negative or
    % both positive.
    s(:) = log(posMin) + prob .* (log(posMax) - log(posMin));
    s(:) = exp(s(:));
    
    if bothNegative
        % Switch back from all positive to all negative.
        s = -s;
    elseif xminNegativeOnly
        % Readjust the result back into the original range, by flipping
        % about posMin.
        s = posMin - (s - posMin);
    end
end

end
