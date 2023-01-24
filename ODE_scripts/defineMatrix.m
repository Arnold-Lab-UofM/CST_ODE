

%%
function matrix = defineMatrix(items, nsamples, itemType)

% items is a cell array. Each element of items is a cell array with 4 elements
% - a name, a distribution string and 2 numbers.
%
% Create a matrix of values for the items.  The matrix will have a column for
% each item and nsamples rows.
%
% For each item without a distribution use its first numeric value as its value
% for all runs, i.e.  for each element of its column.
%
% For each item with a distribution choose a set of random values according to
% the distribution, using its 2 numeric values as distribution parameters (ex.
% min/max, mean/stddev).

itemCount = length(items);
matrix = zeros(nsamples, itemCount);

for i = 1:itemCount
    item = items{i};
    matrix(:,i) = getSampleValues(item, nsamples, itemType, i);
end

end

function sampleValues = getSampleValues(item, nsamples, itemType, itemOrdinal)
%
% Get a vector of nsamples number of values for an item based on its
% distribution code.

if (isempty(item{2}))
    sampleValues = ones(1, nsamples) * item{3};
elseif (strcmp(item{2}, 'u'))
    sampleValues = lhs_ode_unif_new(item{3}, item{4}, nsamples, false);
elseif (strcmp(item{2}, 'lu'))
    sampleValues = lhs_ode_unif_new(item{3}, item{4}, nsamples, true);
elseif (strcmp(item{2}, 'n'))
    sampleValues = lhs_ode_norm_new(item{3}, item{4}, nsamples);
else
    msgIdent = 'LHS:InvalidDistribution';
    msgfmt = '%s %d %s has invalid distribution ''%s''';
    msg = sprintf(msgfmt, itemType, itemOrdinal, item{1}{1}, item{1}{2});
    throw(MException(msgIdent, msg));
end

end