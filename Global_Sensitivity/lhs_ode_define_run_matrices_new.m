function  [paramMatrix, icMatrix, pulseMatrix, LHSmatrix, PRCC_var] = ...
            lhs_ode_define_run_matrices_new(odeSettings)

% Inputs:
%
% odeSettings is a structure containing the settings for an ODE LHS.
%
% odeSettings.NR is the number of rows in the parameters and initialConditions
% matrices.
%
% odeSettings.parameters are the parameters for the model. Each element
% specifies whether or not the parameter is to be varied, and if so what
% probability distribution to use.
%
% odeSettings.initialConditions are the initial conditions for the model.  Each
% element specifies whether or not the initial condition is to be varied, and
% if so what probability distribution to use.
%
% Parameters and initial conditions to be varied are defined here using the LHS
% (Latin Hypercube Sampling) method.
%
% Once varied parameters and initial conditions are defined, any computed
% parameters and computed initial conditions are defined.
%
% odeSettings.pulseTimePoints is a vector of time points when to pulse the
% initial conditions.
%
% odeSettings.pulseValues is a cell array for specifying pulse values. It is in
% a similar format as odeSettings.parameters and odeSettings.initialConditions.
%
% Outputs:
%
% paramMatrix has a column for each parameter in odeSettings.parameters. The
% number of rows is the number of model runs to be performed, odeSettings.NR.
%
% icMatrix has a column for each initial condition in
% odeSettings.initialConditions.  The number of rows is the number of model
% runs to be performed, odeSettings.NR.
%
% pulseMatrix(NR, NTP, NIC) is a matrix of pulse values derived from
% odeSettings.pulseValues. pulseMatrix(r, tp, ic) is the pulse value for
% initial condition ic at pulse time point tp for run r. NR is the number of
% runs, NTP is the number of time points and NIC is the number of initial
% conditions. If an initial condition, ic, is not being pulsed then
% pulseMatrix(:, :, ic) is 0.
%
% LHSmatrix is the columns of paramMatrix for the varied parameters, the
% columns of icMatrix for the varied initial conditions and the columns of
% pulseMatrix for the initial conditions pulsed, for each pulse time point.
%
% PRCC_var is a cell array of the names of the varied parameters, initial
% conditions and pulse items. This is needed by code which performs PRCC
% analysis. PRCC analysis determines which of the varied items are significant.
% That code needs to label the results with names of the significant varied
% items.  It makes for simpler PRCC analysis code if just the list of varied
% item names is provided, rather than requiring the PRCC analysis code to
% perform multiple iterations over the parameter, initial condition and pulse
% cell arrays to find those names.
% 

paramMatrix = defineMatrix(odeSettings.parameters, odeSettings.NR, ...
                           'parameter');

icMatrix = defineMatrix(odeSettings.initialConditions, odeSettings.NR, ...
                        'initial condition');

% Update the parameter and initial condition matrices for any computed
% parameters and any computed initial conditions.
if (isfield(odeSettings, 'computeParamsInitCondHandle'))
    for i = 1:odeSettings.NR
        [paramMatrix(i,:), icMatrix(i,:)] = ...
            odeSettings.computeParamsInitCondHandle(paramMatrix(i,:), ...
                                                    icMatrix(i,:));
    end
end

% Construct an LHS matrix to use for performing PRCC analysis.  It is the
% columns of the parameters matrix for those parameters that are varied,
% and the colums of the initial condition matrix for those initial
% conditions that are varied.
LHSmatrix = [];
PRCC_var = {};
[LHSmatrix, PRCC_var] = defineLHSmatrix(odeSettings.parameters, ...
                                        paramMatrix, LHSmatrix, PRCC_var);

[LHSmatrix, PRCC_var] = defineLHSmatrix(odeSettings.initialConditions, ...
                                        icMatrix, LHSmatrix, PRCC_var);


% Define the pulse matrix. This must come after defining the computed
% initial conditions, since those might be used to define parts of the
% pulse matrix. And it must come after the LHSMatrix and PRCC_var are
% updated based on the varying parameters and initial conditions.
pulseMatrix = [];
if isfield(odeSettings, 'pulseTimePoints')
    pulseCount = length(odeSettings.pulseTimePoints);
    [pulseMatrix, LHSmatrix, PRCC_var] = ...
        definePulseMatrix(pulseCount, odeSettings.pulseValues, ...
        odeSettings.NR, odeSettings.icNames, icMatrix, LHSmatrix, PRCC_var);
end

end

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

function [LHSmatrix, PRCC_var] = defineLHSmatrix(items, itemMatrix, ...
                                                 LHSmatrix, PRCC_var)
% Udpate LHSmatrix and PRCC_var with information for the varied items.
for i = 1:length(items)
    name = items{i}{1};
    distribution = items{i}{2};
    if (~isempty(distribution))
        % Item has a distribution, so it is varied.
        LHSmatrix(:,end+1) = itemMatrix(:,i);
        PRCC_var{end+1} = name;
    end
end

end % function defineLHSMatrix

function [pulseMatrix, LHSmatrix, PRCC_var] = definePulseMatrix( ...
    pulseTimePointCount, pulseValues, nsamples, icNames, icMatrix, ...
    LHSmatrix, PRCC_var)

% Define the pulse matrix. Update the LHSMatrix and PRCC_var based on the
% pulse items that are varied.

pulseMatrix = zeros(nsamples, pulseTimePointCount, length(icNames));
pulseOrdinal = 0;

% How many times each initial condition is referenced in a pulse value.
% At the end it should be exactly pulseTimePointCount for those that are
% referenced in a pulse value, and 0 otherwise.
%
% It is also the ordinal of the time point for the current pulse value
% item. This specifies the time point to update in the pulse matrix for the
% current pulse value.
icRefCount = zeros(length(icNames));

for item = pulseValues
    pulseOrdinal = pulseOrdinal + 1;
    name = item{1}{1};
    idx = getIcIndex(name, icNames, pulseOrdinal);

    [rows, columns] = size(item{1});
    if rows > 1
        msgIdent = 'LHS:ODE:RUN:InvalidSetting';
        msgfmt = 'Cell array for pulse value %d has more than 1 row, %d';
        msg = sprintf(msgfmt, pulseOrdinal, rows);
        throw(MException(msgIdent, msg));     
    end
    
    if columns == 1
       % Use the values of the initial condition for all runs. This means
       % the pulse value is not being varied. Even if the initial condition
       % to be pulsed is being varied, the pulse values are not being
       % varied independently, so should not be included in a PRCC
       % analysis.
       icRefCount(idx) = icRefCount(idx) + 1;
       checkIcRefCount(icRefCount(idx), pulseTimePointCount, pulseOrdinal, ...
                       name);
       pulseMatrix(:,icRefCount(idx),idx) = icMatrix(:,idx);
    elseif columns == 4
        icRefCount(idx) = icRefCount(idx) + 1;
        checkIcRefCount(icRefCount(idx), pulseTimePointCount, pulseOrdinal, ...
                        name);
        sampleValues = getSampleValues(item{1}, nsamples, 'pulse value');
        pulseMatrix(:,icRefCount(idx),idx) = sampleValues;
        
        if ~ isempty(item{1}{2})
            % The pulse value is varied.
            LHSmatrix(:,end+1) = sampleValues;
            
            % Don't reuse the ic name, modify it so it's unique in
            % PRCC_var.
            PRCC_var{end+1} = sprintf('%s_pulse_%d',name, icRefCount(idx));
        end
    else
        msgIdent = 'LHS:ODE:RUN:InvalidSetting';
        msgfmt = ['Cell array for pulse value %d has invalid number of', ...
                 ' elements %d. There should be 1 or 4.'];
        msg = sprintf(msgfmt, pulseOrdinal, columns);
        throw(MException(msgIdent, msg));     
    end   
end

for ic = 1:length(icRefCount)
    
    refCount = icRefCount(ic);
    if refCount ~= 0 && refCount ~= pulseTimePointCount
        msgIdent = 'LHS:ODE:RUN:InvalidIcRefCount';
        msgfmt = ['The number of pulse values for initial condition %s' ...
                  ' was %d. It should be either 0 or %d.'];
        msg = sprintf(msgfmt, char(icNames(ic)), refCount, pulseTimePointCount);
        throw(MException(msgIdent, msg));     
    end
    
end

end

function idx = getIcIndex(name, icNames, pulseOrdinal)

idx = find(ismember(icNames, name));
if isempty(idx)
    msgIdent = 'LHS:ODE:RUN:InvalidPulseICName';
    msgfmt = ['Cell array for pulse value %d has invalid initial condition' ...
              ' name %s'];
    msg = sprintf(msgfmt, pulseOrdinal, name);
    throw(MException(msgIdent, msg));
end

end % function checkIcName

function checkIcRefCount(icRefCount, pulseTimePointCount, pulseOrdinal, icName)

if icRefCount > pulseTimePointCount
    msgIdent = 'LHS:ODE:RUN:InvalidIcRefCount';
    msgfmt = ['Cell array for pulse value %d is the %d reference to' ... 
              ' initial condition %s, which is more than the number of', ...
              ' pulse time points, %d'];
    msg = sprintf(msgfmt, pulseOrdinal, icRefCount, icName, ...
                  pulseTimePointCount);
    throw(MException(msgIdent, msg));
end

end

