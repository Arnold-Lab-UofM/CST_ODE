function ode = lhs_ode_get_run_settings_new(lhsSettingsFileName)

% Get the ODE LHS settings from a settings file.
%
% The caller will perform checks on some of the settings, since
% different callers will have different settings requirements in some
% cases.

if isempty(lhsSettingsFileName)
    title = 'Select ODE LHS settings file';
    default = 'lhs_ode_settings.m';
    [lhsSettingsFileName, lhsDirectory, ~] = uigetfile('*.m', title, default);
    
    % filesep is a builtin Matlab function that returns the file system file
    % separator for the local system. For example, on a Unix system it returns
    % '/'.
    lhsSettingsPathName = [lhsDirectory, filesep, lhsSettingsFileName];

else
    % The directory is not a separate item if specified as a function
    % argument.
    lhsDirectory = '';
    lhsSettingsPathName = lhsSettingsFileName;
end    

if lhsSettingsFileName == 0
    msgIdent = 'LHS:NoSettingsFile';
    msg = 'No ODE LHS settings file specified.\n';
    throw(MException(msgIdent, msg));
end

% Load the LHS settings.
try
    if ~isempty(lhsDirectory)
        curDir = cd(lhsDirectory);
    end
    [~, name, ~] = fileparts(lhsSettingsFileName);
    
    % try to feed thru handle
    settingsHandle = str2func(name);
    ode = settingsHandle();

    if ~isempty(lhsDirectory)
        cd(curDir);
    end
catch e
    msgIdent = 'LHS:ODE:SettingsFileEvalError';
    msg = ['Unable to load ODE LHS settings from file "%s".\n%s', ...
            '\nThis might be due to a typo or', ...
            ' stray characters in the settings file.', ...
            '\nThe settings file must contain valid Matlab code.'];
    throw(MException(msgIdent, msg, lhsSettingsPathName, e.message));
end

assignin('base','lhsOdeSettingsFileName', lhsSettingsFileName);
assignin('base', 'lhsDirectory', lhsDirectory);

% Perform some consistency checks on the loaded settings.

if (~exist('ode', 'var'))
    msgIdent = 'LHS:ODE:MissingSetting';
    msg = 'The ode parameter structure does not exist.';
    throw(MException(msgIdent, msg));    
end

isFieldExist(ode, 'NR');
is_numeric_min(ode.NR, 1, 'NR', 'LHS:ODE:RUN:InvalidSetting');

% Retrieve tspan from workspace
isFieldExist(ode, 'tspan');
if (~isnumeric(ode.tspan))
    msgIdent = 'LHS:ODE:LHS:InvalidTspanList';
    msg = 'The tspan list of model time points is not numeric.';
    throw(MException(msgIdent, msg));  
end

% tspan must be a montonically increasing sequence of numeric values.
% The values in tspan will have a Matlab type of double.
t_previous = -1;
for t = ode.tspan
    if t < 0
        msgIdent = 'LHS:ODE:InvalidTimePoint';
        msgfmt = 'The tspan time point value of %d is < 0';
        msg = sprintf(msgfmt, t);
        throw(MException(msgIdent, msg));       
    elseif t <= t_previous
        msgIdent = 'LHS:ODE:InvalidTimePoint';
        msgfmt = ['The tspan time point value of %d' ...
                    ' is <= the previous time point of %d'];
        msg = sprintf(msgfmt, t, t_previous);
        throw(MException(msgIdent, msg));                   
    end
    
    t_previous = t;
end

% analysisTimePoints must be a montonically increasing sequence of numeric
% values.  The values in analysisTimePoints will have a Matlab type of double.
%
% Each element of time points must also be an element of tspan.
%
% Also define the indices into tspan corrsponding to each analysisTimePoints
% element. These will be used to index into result matrices which will have one
% row for each element in tspan.
%
% For example, suppose tspan = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0], and
% analysisTimePoints = [1.5, 2.5] then time_index = [3, 5], since 3 is the
% index in tspan for 1.5, and 5 is the index in tspan for 2.5. Result matrices
% would have 6 rows, one for each element of tspan. If we want to perform
% analysis on the rows corresponding to time points 1.5 and 2.5, that
% corresponds to using rows 3 and 5 of those result matrices.

isFieldExist(ode, 'analysisTimePoints');

t_previous = -1;
ode.time_index = zeros(1, length(ode.analysisTimePoints));
i = 1; % The next element of ode.time_index to be defined.
for t = ode.analysisTimePoints
    if t < 0
        msgIdent = 'LHS:ODE:InvalidTimePoint';
        msgfmt = 'The analysisTimePoints time point value of %d is < 1';
        msg = sprintf(msgfmt, t);
        throw(MException(msgIdent, msg));       
    elseif t <= t_previous
        msgIdent = 'LHS:ODE:InvalidTimePoint';
        msgfmt = ['The analysisTimePoints time point value of %d' ...
                    ' is <= the previous time point of %d'];
        msg = sprintf(msgfmt, t, t_previous);
        throw(MException(msgIdent, msg));                   
    end
    
    [flag, tspan_index] = ismember(t, ode.tspan);
    if flag == 1
        ode.time_index(i) = tspan_index;
    else
        msgIdent = 'LHS:ODE:InvalidTimePoint';
        msgfmt = 'The analysisTimePoints time point value of %d is not in tspan';
        msg = sprintf(msgfmt, t);
        throw(MException(msgIdent, msg));
    end
    t_previous = t;
    i = i + 1;
end

if (isfield(ode, 'alpha'))
    if (~isnumeric(ode.alpha))
        msgIdent = 'LHS:ODE:InvalidAlpha';
        msgfmt = 'The alpha value of ''%s'' is not numeric.'
        msg = sprintf(msgfmt, ode.alpha);
        throw(MException(msgIdent, msg));
    end
    
    if (ode.alpha <= 0)
        msgIdent = 'LHS:ODE:InvalidAlpha';
        msgfmt = 'The alpha value of %d is <= 0.'
        msg = sprintf(msgfmt, ode.alpha);
        throw(MException(msgIdent, msg));
    end
end

isFieldExist(ode, 'odeModelHandle');
check_function_handle(ode.odeModelHandle, 'odeModelHandle', ...
                      'LHS:ODE:MissingSetting');

if (isfield(ode, 'computeParamsInitCondHandle'))
    check_function_handle(ode.computeParamsInitCondHandle, ...
                         'computeParamsInitCondHandle', ...
                          'LHS:ODE:MissingSetting');
end

isFieldExist(ode, 'solverTimeLimit');
if ode.solverTimeLimit <= 0.0
    msgIdent = 'LHS:ODE:SETTINGS:BadSolverTimeLimit';
    msg = 'solverTimeLimit of %f is <= 0.0.\n';
    error(msgIdent, msg, ode.solverTimeLimit);
end

isFieldExist(ode, 'saveOnError');
if ode.saveOnError ~= 0 && ode.saveOnError ~= 1
    msgIdent = 'LHS:ODE:SETTINGS:BadSaveOnError';
    msg = 'saveOnError of %d is not 0 or 1.\n';
    error(msgIdent, msg, ode.saveOnError);
end

% Check the parameters setting.
isFieldExist(ode, 'parameters');
[varyingParameterCount, ode.paramNames] = checkItems(ode.parameters, ...
                                                     'parameter');

% Check the initial conditions setting.
isFieldExist(ode, 'initialConditions');
[varyingIcCount, ode.icNames] = checkItems(ode.initialConditions, ...
                                          'initial condition');
                                      
ode.varyingItemCount = varyingParameterCount + varyingIcCount;                                     

% Check the output labels.
isFieldExist(ode, 'outputLabels');
icCount = length(ode.initialConditions);
olCount = length(ode.outputLabels);
if (olCount ~= icCount)
    msgIdent = 'LHS:ODE:InvalidSetting';
    msgfmt = ['The number of output labels, %d, does not match the number' ...
             ' of initial conditionsi %d']
    msg = sprintf(msgfmt, olCount, icCount);
    throw(MException(msgIdent, msg));
end


% Check the pulse items, if defined.
% All pulse items should be defined or none should be defined.

definedPulseSettingsCount = 0;
definedPulseSettings = '';
if (isfield(ode, 'pulseTimePoints'))
    definedPulseSettingsCount = definedPulseSettingsCount + 1;
    definedPulseSettings = [definedPulseSettings 'pulseTimePoints'];
end

if (isfield(ode, 'pulseValues'))
    definedPulseSettingsCount = definedPulseSettingsCount + 1;
    definedPulseSettings = [definedPulseSettings 'pulseValues'];
end

if (definedPulseSettingsCount ~= 0 && definedPulseSettingsCount ~= 2)
    msgIdent = 'LHS:ODE:InvalidPulseSettings';
    msgfmt = 'Some but not all of the pulse settings are defined: %s.';
    msg = sprintf(msgfmt, definedPulseSettings);
    throw(MException(msgIdent, msg));       
end

if (definedPulseSettingsCount == 2)

    % pulseTimePoints should not be empty.
    if isempty(ode.pulseTimePoints)
        msgIdent = 'LHS:ODE:InvalidPulseTimesteps';
        msg = 'The pulseTimePoints vector is empty';
        throw(MException(msgIdent, msg));       
    end

    % pulseTimePoints must be a montonically increasing sequence of numeric
    % values. The values will have a Matlab type of double. 
    % Also each pulse time point should be in the tspan vector.
    pts_previous = -1;
    ode.pulseTimePointIdx = [];
    for pts = ode.pulseTimePoints
        if pts < 0
            msgIdent = 'LHS:ODE:InvalidPulseTimePoint';
            msgfmt = 'The pulse time point value of %d is < 0';
            msg = sprintf(msgfmt, pts);
            throw(MException(msgIdent, msg));       
        elseif pts <= pts_previous
            msgIdent = 'LHS:ODE:InvalidPulseTimePoint';
            msgfmt = ['The pulse time point value of %d' ...
                        ' is <= the previous time point of %d'];
            msg = sprintf(msgfmt, pts, pts_previous);
            throw(MException(msgIdent, msg));                   
        end

        pulseTimePointIdx = find(ode.tspan == pts);
        if isempty(pulseTimePointIdx)
            msgIdent = 'LHS:ODE:SETTINGS:BadPulseTimestep';
            msg = 'pulseTimestep of  %d not found in tspan\n';
            error(msgIdent, msg, pts);
        else
            if pulseTimePointIdx == 1
                msgIdent = 'LHS:ODE:SETTINGS:BadPulseimePointIdx';
                msg = 'pulseTimePointIdx is 1. Can''t pulse on the 1st time step.\n';
                error(msgIdent, msg);
            end
        end
        
        ode.pulseTimePointIdx = [ode.pulseTimePointIdx pulseTimePointIdx];
        
        pts_previous = pts;
    end

    if (isempty(ode.pulseValues))
        msgIdent = 'LHS:ODE:InvalidPulseValues';
        msg = 'The pulseValues vector is empty';
        throw(MException(msgIdent, msg));       
    else
        ode.varyingItemCount = checkPulseValues(ode.pulseValues, ...
                                    ode.icNames, ode.varyingItemCount);
    end

end % if (definedPulseSettingsCount == 2)

if ode.varyingItemCount == 0
    % Not varying any items. Should only do 1 run in this case.
    if (ode.NR ~= 1)
        msgIdent = 'LHS:ODE:RUN:InvalidSetting';
        msgfmt = ['There are no varying items (parameters or initial' ...
                  ' conditions) but the number of runs (NR) of %d is not' ...
                  ' equal to 1. It does not make sense to do more than one' ...
                  ' run if not varying any items.'];
        msg = sprintf(msgfmt, ode.NR);
        throw(MException(msgIdent, msg));                   
    end
elseif ode.varyingItemCount == 1
    % Not useful to vary only 1 parameter - treat as an error.
    msgIdent = 'LHS:ODE:InvalidSetting';
    msg = ['There is only 1 varying imte (parameter or initial condition.' ...
           ' PRCC cannot be performed on only 1 varying item (i.e. for 2' ... 
           ' or more runs) and it does not make sense to vary 1 parameter' ...
           ' for only 1 run.'];
    throw(MException(msgIdent, msg));
else
    % Varying more than 1 parameter, must do at least 2 runs.
    is_numeric_min(ode.NR, 2, 'NR', 'LHS:ODE:RUN:InvalidSetting');
end


end

function isFieldExist(ode, field)

if (~isfield(ode, field))
    msgIdent = 'LHS:ODE:MissingSetting';
    msgfmt = '%s setting does not exist.';
    msg = sprintf(msgfmt, field);
    throw(MException(msgIdent, msg));    
end

end

function [varyingCount, itemNames] = checkItems(items, itemTypeName)
%
% Items is a cell array.
%
% Each element of items is a cell array that should match the format checked by
% function checkItem.
%
% varyingCount: The number of items that are being varied.
%
% itemNames: A list of just the names of the items. Useful when some other
%            data refers to these items by name, to check that the
%            reference is to a valid item name.

varyingCount = 0;
itemCount = 0;
itemNames = {};
for item = items
    itemCount = itemCount + 1;
    [varyingCount, itemNames] = checkItem(item, itemTypeName, itemCount, ...
                                          varyingCount, itemNames);
end % end for

end % function checkItems

function [varyingCount, itemNames] = checkItem(item, itemTypeName, ...
                                        itemCount, varyingCount, itemNames)
%
% Check a single item.
%
% Each item should have the following format.
%
% A string for the item name.
%
% A string for a probability distribution.
%     '':   no distribution.
%     'u':  uniform distribution
%     'lu': log uniform distribution
%     'n':  normal distribution
%
% 2 numeric values.
%
% If the distribution is '' then the 2 numeric values should be the same.
%
% If the distribution is 'u' or 'lu' then the 2 numeric values are the min and
% max to use for the distribution and should be different and min < max.
%
% If the distribution is 'n' then the 2 numeric values are the mean and stddev
% to use for the distribution and stddev > 0.0.
%
% name: A string used in error messages. It is for the type of items being
%       processed, ex. 'parameter' or 'initial condition'.
%
% variedCount: A count of the varied items in a list of items. It is
%              incremented if this item is varied (has a valid distribution).
%
% itemNames: A list pf names for the items. It has this item's name appended.
%
% The return values are:
%
% varyingCount: Incremented if this items is one that is being varied.
%
% itemNames: Updated with this item's name.
%

% A version of name with the first character in upper case, when the name is to
% be used at the start of an error message.
upperItemTypeName = [upper(itemTypeName(1)), itemTypeName(2:end)];

% The item must be exactly 4 items, 2 strings and 2 numbers.

[rows, columns] = size(item{1});
if (rows > 1)
    msgIdent = 'LHS:ODE:RUN:InvalidSetting';
    msgfmt = 'Cell array for %s %d has more than 1 row, %d';
    msg = sprintf(msgfmt, itemTypeName, itemCount, rows);
    throw(MException(msgIdent, msg));
    
end

if (columns ~= 4)
    msgIdent = 'LHS:ODE:RUN:InvalidSetting';
    msgfmt = 'Cell array for %s %d has %d columns, instead of 4';
    msg = sprintf(msgfmt, itemTypeName, itemCount, columns);
    throw(MException(msgIdent, msg));
end

if (~ischar(item{1}{1}))
    msgIdent = 'LHS:ODE:RUN:InvalidSetting';
    msgfmt = ['1st element of cell array for %s %d' ...
        ' is not a character string'];
    msg = sprintf(msgfmt, itemTypeName, itemCount);
    throw(MException(msgIdent, msg));
end

if (~ischar(item{1}{2}))
    msgIdent = 'LHS:ODE:RUN:InvalidSetting';
    msgfmt = ['2nd element of cell array for %s %d %s' ...
        ' is not a character string'];
    msg = sprintf(msgfmt, itemTypeName, itemCount, item{1}{1});
    throw(MException(msgIdent, msg));
end

if (~isnumeric(item{1}{3}))
    msgIdent = 'LHS:ODE:RUN:InvalidSetting';
    msgfmt = ['3rd element of cell array for %s %d %s' ...
        ' is not numeric'];
    msg = sprintf(msgfmt, itemTypeName, itemCount, item{1}{1});% Pulse items don't always have exactly the same format as parameters or initial
% conditions. They might have 4 items per entry like they do, or only 1.
    throw(MException(msgIdent, msg));
end

if (~isnumeric(item{1}{4}))
    msgIdent = 'LHS:ODE:RUN:InvalidSetting';
    msgfmt = ['4th element of cell array for %s %d %s' ...
        ' is not numeric'];
    msg = sprintf(msgfmt, itemTypeName, itemCount, item{1}{1});
    throw(MException(msgIdent, msg));
end

itemNames{end+1} = item{1}{1};

distribution = item{1}{2};

if (isempty(distribution))
    if (item{1}{3} ~= item{1}{4})
        msgIdent = 'LHS:ODE:RUN:InvalidSetting';
        msgfmt = ['%s %s does not have a distribution,' ...
            ' but the 2 numeric values are not the same: %f %f'];
        msg = sprintf(msgfmt, upperItemTypeName, item{1}{1}, item{1}{3}, item{1}{4});
        throw(MException(msgIdent, msg));
    end
elseif (strcmp(distribution, 'u') || strcmp(distribution, 'lu'))
    if (item{1}{3} == item{1}{4})
        msgIdent = 'LHS:ODE:RUN:InvalidSetting';
        msgfmt = ['%s %s has a distribution of ''%s'',' ...
            ' but the min and max are the same: %f %f'];
        msg = sprintf(msgfmt, upperItemTypeName, item{1}{1}, item{1}{2}, ...
            item{1}{3}, item{1}{4});
        throw(MException(msgIdent, msg));
    elseif (item{1}{3} > item{1}{4})
        msgIdent = 'LHS:ODE:RUN:InvalidSetting';
        msgfmt = ['%s %s has a distribution of ''%s'',' ...
            ' but min, %f, >  max, %f'];
        msg = sprintf(msgfmt, upperItemTypeName, item{1}{1}, item{1}{2}, ...
            item{1}{3}, item{1}{4});
        throw(MException(msgIdent, msg));
    end
    varyingCount = varyingCount + 1;
    
elseif(strcmp(distribution, 'n'))
    if (item{1}{4} <= 0.0)
        msgIdent = 'LHS:ODE:RUN:InvalidSetting';
        msgfmt = ['%s %s has a distribution of ''n'',' ...
            ' but the stddev of %f is <= 0.0.'];
        msg = sprintf(msgfmt, upperItemTypeName, item{1}{1}, item{1}{4});
        throw(MException(msgIdent, msg));
    end
    varyingCount = varyingCount + 1;
else
    msgIdent = 'LHS:ODE:RUN:InvalidSetting';
    msgfmt = ['%s %s has invalid distribution ''%s''.' ...
        ' It must be ''u'', ''lu'' or ''n'''];
    msg = sprintf(msgfmt, upperItemTypeName, item{1}{1}, item{1}{2});
    throw(MException(msgIdent, msg));
end
end


function varyingCount = checkPulseValues(pulseValues, icNames, varyingCount)
% Check the pulse values for validity.
%
% pulseValues is a cell array of pulse value items.
%
% Each pulse value item is a cell array with either 4 elements, in the same
% format as a parameter or initial condition, or with one element. In
% either case the 1st element must be the same is the name of an initial
% condition - the initial condition its values are to pulse.
%
% If only the initial condition name is specified then the pulse is the
% same as the value of the initial condition.
%
% If 4 items are specified, then they define the value of the pulse. If the
% 2nd item is empty (no distribution specified) then a single value for is
% used for the pulse all runs. If the 2nd item is not empty (a valid
% distribution is specified) then a set of values is used, with a different
% value for each run.
% 
% Each initial condition to be pulsed must have as many pulse value items as
% there are pulse time points. This allows using a different pulse value,
% or set of pulse values, for each pulse time point.
%return;
itemCount = 0;

% Argument required by function checkItem, but not needed here.
dummyNames = {};

for item = pulseValues
    itemCount = itemCount + 1;

    [rows, columns] = size(item{1});
    if rows > 1
        msgIdent = 'LHS:ODE:RUN:InvalidSetting';
        msgfmt = 'Cell array for pulse value %d has more than 1 row, %d';
        msg = sprintf(msgfmt, itemCount, rows);
        throw(MException(msgIdent, msg));     
    end
    
    if columns == 1
        % Check ic name.
        checkIcName(item{1}{1}, icNames, itemCount);
    elseif columns == 4
        % Check ic name and then check all the elements of the item.
        checkIcName(item{1}{1}, icNames, itemCount);
        [varyingCount, ~] = checkItem(item, 'pulseValue', itemCount, ...
                                      varyingCount,  dummyNames);
    else
        msgIdent = 'LHS:ODE:RUN:InvalidSetting';
        msgfmt = 'Cell array for pulse value %d has has invalid number of'
                 ' elements %d. There should be 1 or 4.';
        msg = sprintf(msgfmt, itemCount, columns);
        throw(MException(msgIdent, msg));     
    end   
end % end for

end % function checkPulseValues

function checkIcName(name, icNames, pulseOrdinal)

idx = find(ismember(icNames, name));
if isempty(idx)
    msgIdent = 'LHS:ODE:RUN:InvalidPulseICName';
    msgfmt = ['Cell array for pulse value %d has invalid initial condition' ...
              ' name %s'];
    msg = sprintf(msgfmt, pulseOrdinal, name);
    throw(MException(msgIdent, msg));
end

end % function checkIcName
