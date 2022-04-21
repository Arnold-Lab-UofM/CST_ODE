% Kirschner Lab, University of Michigan
% Kirschner Lab website: http://malthus.micro.med.umich.edu/lab/

% Last edit: March 11th, 2020 by
% Arnold Lab, University of Michigan
% Melissa Lemke, PhD Candidate
% Modified for Ab-FcR model - changed the ODE solver to ode15s lines 135, 149, 187

function lhs_ode_run_new(lhsODESettingsFileName)

% -------------------------------------------------------------------------
%
% Setup and run of an LHS analysis of an ODE model - main routine
%
% The results should be compared to the PRCC results section in
% Supplementary Material D and Table D.1 for different N (specified by
% "NR" in the script below
%
% Input:
%
% lhsOdeSettingsFileName: a file that contains the definitions of the LHS
%                         settings needed by this script. If not specified
%                         as a function argument then retrieved using a
%                         file browser dialog.
%
% lhsDirectory: The main LHS directory to contain LHS settings file.
%                 Defined by the selection made in the Matlab file browser
%                 dialog when selecting the LHS settings file.
%
% See the header comments in a settings file for a description of the members
% of the returned structure named "odeSettings".
%
%
% Output
%
% The LHS analysis results are displayed in the Matlab workspace and can be
% saved to a binary file. PRCC analyses can also be run on the LHS analysis
% results vectors.
%
%--------------------------------------------------------------------------
% 
% Declare variables that can be saved to a .mat file. They must be also
% assigned to the base workspace, otherwise the cleanup function that does
% the saving won't be able to access them. Each variable to saved must be
% declared and assigned here, since not all will be defined later,
% depending on program options (ex. a PRCC might not be performed). Also
% each variable must be assigned again if it is defined later, so the most
% current value gets saved - assigin essentially copies the current value
% to the base workspace, it doesn't cause any future values to be available
% in the base workspace.
odeSettings = {}; assignin('base', 'odeSettings', odeSettings);
odeOut = {}; assignin('base', 'odeOut', odeOut);
totalPrcc = []; assignin('base', 'totalPrcc', totalPrcc);
totalSign = []; assignin('base', 'totalSign', totalSign);
totalPrccSignEval = []; assignin('base', 'totalPrccSignEval', ...
                                 totalPrccSignEval);
paramMatrix = []; assignin('base', 'paramMatrix', paramMatrix);
icMatrix = []; assignin('base', 'icMatrix', icMatrix);
pulseMatrix = []; assignin('base', 'pulseMatrix', pulseMatrix);
LHSmatrix = []; assignin('base', 'LHSmatrix', LHSmatrix);
PRCC_var = []; assignin('base', 'PRCC_var', PRCC_var);

                                
% This is for handling the case where the Matlab ode solver fails for a
% set of parameters and initial conditions. A solver can fail if it can't
% find a sufficiently small time step, if it encounters a singular matrix
% or if it takes longer than our defined solver time limit.
error_runs={};
assignin('base', 'error_runs', error_runs);                                 

cleanupObj = onCleanup(@cleanMeUp);

fprintf('\n=====================================\n');
fprintf('   ODE LHS Run');
fprintf('\n=====================================\n');


% Load the settings file - everything returned in the odeSettings structure.
odeSettings = lhs_ode_get_run_settings_new(lhsODESettingsFileName);

%Melissa Lemke edit - add extra parameter sets in case of error
odeSettings.NR_original = odeSettings.NR;
odeSettings.NR = odeSettings.NR_original+ceil(odeSettings.NR_original*.1);

% Put current value in base workspace so the cleanup function can save it.
assignin('base', 'odeSettings', odeSettings);

% Build the parameter matrix and initial condition matrix.
%
% The ith row of the paramter matrix are the parameters to use for run i.
%
% The ith row of the initial condition matrix are the initial conditions to use
% for run i.

[paramMatrix, icMatrix, pulseMatrix, LHSmatrix, PRCC_var] = ...
    lhs_ode_define_run_matrices_new(odeSettings);

assignin('base', 'paramMatrix', paramMatrix);
assignin('base', 'icMatrix', icMatrix);
assignin('base', 'pulseMatrix', pulseMatrix);
assignin('base', 'LHSmatrix', LHSmatrix);
assignin('base', 'PRCC_var', PRCC_var);

% The function to call after each successful solver step. Needed to enforce
% the ODE solver time limit.
outputFun = @(t,y,flag)solverTimeLimitFun(t, y, flag, odeSettings.solverTimeLimit);
options = odeset('OutputFcn', outputFun,'AbsTol',1e-50,'RelTol',1e-10);%

% Pre-allocate arrays for speed based on the output variable list 
odeOut = initVars(length(odeSettings.initialConditions), ...
                  length(odeSettings.tspan), odeSettings.NR_original);

% These will be used to contain warnings from ODE solver failures, including
% our ODE solver time limit facility.
savedWarnMsg = '';
savedWarnId = '';

% These might not be defined if an ODE solver error occurs before the first
% successful integration step. Our error handling requires that they be
% defined.
t = [];
y = [];

% Perform the NR model runs.
% for run = 1:odeSettings.NR - Melissa Lemke comment out
%Melissa Lemke edit
run = 1;
sucessful_run = zeros(1,odeSettings.NR);

while sum(sucessful_run) < odeSettings.NR_original %end of Melissa Lemke edit: loop starts at zero, so "less than" is actually "less than or eq to"
    fprintf('%d ', run);
    
    if (mod(run, 25) == 0)
        fprintf('\n');
    end

    if run>odeSettings.NR
        fprintf('\n');
        fprintf('Not enough successful runs');
        fprintf('\n');
        return
    end
    
    runParameters = paramMatrix(run,:);
    y0 = icMatrix(run,:);
    
    
    try
        if ~isfield(odeSettings, 'pulseTimePoints')
            % Not pulsing ODE conditions, just run the model.
            lastwarn('') % Clear last warning message
            [t,y] = ode15s(@(t,y)odeSettings.odeModelHandle(t, y, runParameters), odeSettings.tspan, y0, options);
            [savedWarnMsg, savedWarnId] = checkWarning();
            
        else
            % Pulsing ODE conditions. Update the initial conditions at the
            % pulse time step.
            
            % For getting pulse amounts. y0 is updated after each solver
            % call, so we need to save the initial conditions, since the
            % pulse amounts may be taken from there.
            y0Initial = y0;
            
            % Solve for the time steps up to just before the first pulse.
            preTspan = odeSettings.tspan(1:odeSettings.pulseTimePointIdx(1) - 1);
            lastwarn('') % Clear last warning message
            [t,y] = ode15s(@(t,y)odeSettings.odeModelHandle(t, y, runParameters), preTspan, y0, options);
            [savedWarnMsg, savedWarnId] = checkWarning();
            
            
            % Only do the pulses if the initial ODE solver call succeeded.
            if (length(t) == length(preTspan))

                % Initial conditions for first pulse solver call, the ODE state
                % at the end of this pre-pulse solver call.
                y0 = y(end,:);

                for i = 1:length(odeSettings.pulseTimePointIdx)
                    
                    % ptpidx and ptpidxNext are indices into tspan.
                    ptpidx = odeSettings.pulseTimePointIdx(i);
                    if i < length(odeSettings.pulseTimePointIdx);
                        % About to do the next pulse. Solve to the time point
                        % just before the time point of the next pulse.
                        ptpidxNext = odeSettings.pulseTimePointIdx(i+1) - 1;
                    else
                        % About to do the last pulse. Solve to the end of the
                        % time points.
                        ptpidxNext = length(odeSettings.tspan);
                    end
                    
                    % Pulse.
                    % pulseMatrix is 3D, so is pulseMatrix(run,i,:), y0 2D.
                    pulseVector = reshape(pulseMatrix(run,i,:), size(y0));
                    y0 = y0 + pulseVector;
                   
                    % Solve the ODE model from the end of the last solution to
                    % the time of the next pulse. Note that the 1st time point
                    % is the same as the last time point for the last solver
                    % call, since the model is not solved for that 1st time
                    % point - the initial conditions passed to the solver are
                    % also the 1st model output, for that 1st time point.
                    currTspan = odeSettings.tspan(ptpidx - 1:ptpidxNext);
                    
                    lastwarn('') % Clear last warning message
                    [tCurr,yCurr] = ode15s(@(tCurr,yCurr)odeSettings.odeModelHandle(tCurr, yCurr, runParameters), currTspan, y0, options);
                    [savedWarnMsg, savedWarnId] = checkWarning();

    
                    % Append the current result vectors into final result
                    % vectors. Skip the first time point and first model
                    % output. These are the same as for the end of the prior
                    % solver call.
                    t = [t;tCurr(2:end)];
                    y = [y;yCurr(2:end,:)];

                    % Check if the ode solver stopped early.
                    if (length(tCurr) ~= length(currTspan))
                        break;
                    end                    
    
                    % The initial conditions for the next solver call are the
                    % ODE state at the last time step (i.e. end) of the prior
                    % solver call.                
                    y0 = yCurr(end,:);
                end % end for i = 1:length(odeSettings.pulseTimePointIdx)
            end
        end
    catch err
        % Save the error info so it can be included in error_runs later.
        savedWarnMsg = err.message;
        savedWarnId = err.identifier;
    end
    
    if (length(t) == length(odeSettings.tspan))
         %Melissa Lemke edit
            sucessful_run(run) = 1; 
            run = run+1;% end of Melissa Lemke edit
            
        % Store ode outputs in the odeOut cell array. Melissa Lemke edit
        odeOut = storeModelRunOutput(odeOut, sum(sucessful_run), y);
           
        % Put current value in base workspace so the cleanup function can save it.
        assignin('base', 'odeOut', odeOut);
    else
        % The ode solver stopped early.
        error_runs{end+1} = {run, savedWarnMsg, savedWarnId};
        assignin('base', 'error_runs', error_runs);
        %Melissa Lemke edit - skip to next run and add to NR if error
        sucessful_run(run) = 0;
        run = run + 1; %end
        
        if (odeSettings.saveOnError == 1) 
            % Pad with last value and record
            odeOut = saveError(odeOut, run, y, odeSettings.tspan);
            
            % Put current value in base workspace so the cleanup function can save it.
            assignin('base', 'odeOut', odeOut);
        end    
    end
end
fprintf('\n');

% Put current value in base workspace so the cleanup function can save it.
assignin('base', 'odeOut', odeOut);

% Do the PRCC, if requested (alpha exists in the ode settings).
%
% Only if doing at least 2 runs and have at least 2 varying parameters.
% If fewer than 2 varying parameters then a fatal Matlab error would be
% generated in the PRCC analysis. Also don't do PRCC if any parameter sets
% had ODE solver errors.

if (~isfield(odeSettings, 'alpha'))
    return;
end

if odeSettings.NR < 2
    msg = ['Skipping PRCC - not enough runs.', ...
           'There are only %d runs but at least 2 are required.\n'
            ];
    fprintf(msg, odeSettings.NR);
    return;
elseif length(PRCC_var) < 2
    msg = 'Skipping PRCC - not enough varying items to do a PRCC.\n';
    fprintf(msg);    
    return;
    %Melissa Lemke comment out and then edit line
% elseif ~ isempty(error_runs)
elseif size(error_runs,1) > odeSettings.NR-odeSettings.NR_original
    msg = 'Skipping PRCC - ODE solver errors occurred for %d parameter sets.\n';
    fprintf(msg, length(error_runs));    
    return;
end

% totalPrcc, totalSign are each 3D matrices: parameters X
% time X outputs.

% That is totalPrcc[1,:,:] is for parameter 1, and is a 2D
% matrix that is time X outputs - 1 row for each time point and 1
% column for each output. 

varyingItemCount = odeSettings.varyingItemCount;
timePointCount = length(odeSettings.analysisTimePoints);
outputCount = length(odeSettings.initialConditions);

totalPrcc = zeros(timePointCount, varyingItemCount, outputCount);
totalSign = zeros(timePointCount, varyingItemCount, outputCount);

eqnCount = length(odeSettings.initialConditions);
fprintf('Calculating PRCC for equation (of %d):\n', eqnCount);
for i = 1:eqnCount

    fprintf(' %d', i);
    if (mod(i, 25) == 0)
        fprintf('\n');
    end


    modelOutput = odeOut{i};

    % Run the PRCC analysis for each output vector in odeOut.
    [prcc prcc_significance] = lhs_ode_prcc_new(LHSmatrix(find(sucessful_run),:), modelOutput, ...
                                                odeSettings.time_index);
    
    % Store the PRCC results in the 3D matrices.
    totalPrcc(:,:,i) = prcc;
    totalSign(:,:,i) = prcc_significance;

end
fprintf('\n');

% Put current values in base workspace so the cleanup function can
% save it.
assignin('base', 'totalPrcc', totalPrcc);
assignin('base', 'totalSign', totalSign);
          
% Define a copy of totalPrcc where each element in the copy is
% set to NAN if the significance for that PRCC value is not
% within the threshold.
totalPrccSignEval = prcc_significance_eval(totalPrcc, totalSign, ...
                                           odeSettings.alpha);

% Put current value in base workspace so the cleanup function can
% save it.
assignin('base', 'totalPrccSignEval', totalPrccSignEval);
        
fprintf('\n=====================================\n');
fprintf('   End ODE LHS Run');
fprintf('\n=====================================\n');

end

function status = solverTimeLimitFun(t,y,flag, solverTimeLimit)   %#ok<INUSL>

% This function gets called after each successful ODE solver integration step.
% It checks the elapsed time from the call to the solver to see if it is over
% the limit.
% Don't throw an error here because that will also trigger the clean up
% function, which will stop processing subsequent parameter sets. We just
% want to halt the solver on those parameter sets that take too long, and
% continue with all the others.

% Instead, issue a warning, since that will be handled by the same error
% handling code that handles warnings from the ODE solver.

persistent INIT_TIME;
status = 0;
switch(flag)
    case 'init'
        INIT_TIME = tic;
    case 'done'
        clear INIT_TIME;
    otherwise
        elapsed_time = toc(INIT_TIME);
        if elapsed_time > solverTimeLimit
            clear INIT_TIME;
            status = 1;
            msg = sprintf('Exceeded time limit of %f, elapsed time is %f,\n',...
                    solverTimeLimit, elapsed_time);
            warning('LHS_ODE_RUN:SolverTimeLimitElapsed', msg);
        end
end
end

function [savedWarnMsg, savedWarnId] = checkWarning()

% Save any non-empty warning so the caller can associate it with any ODE solver
% error.

savedWarnMsg = '';
savedWarnId = '';
[warnMsg, warnId] = lastwarn;
if ~isempty(warnMsg)
    savedWarnMsg = warnMsg;
    savedWarnId = warnId;
end            
end

function cleanMeUp()

% The items from the main function (lhs_ode_run_new) that are to be saved.
% The main function must put these in the base workspace first, otherwise
% they will not be accessible here. Specifically, they are not in the
% caller workspace, since the main function is not the direct caller of
% this clean up function.

odeSettings = evalin('base', 'odeSettings');
odeOut = evalin('base', 'odeOut');
totalPrcc = evalin('base', 'totalPrcc');
totalSign = evalin('base', 'totalSign');
totalPrccSignEval = evalin('base', 'totalPrccSignEval');

paramMatrix = evalin('base', 'paramMatrix');
icMatrix = evalin('base', 'icMatrix');
pulseMatrix = evalin('base', 'pulseMatrix');
LHSmatrix = evalin('base', 'LHSmatrix');
PRCC_var = evalin('base', 'PRCC_var');

error_runs = evalin('base', 'error_runs');

[file, path] = uiputfile('*.mat', 'Save Workspace As', 'Model_LHS.mat');

% Make sure a file name was returned.
if (ischar(file))
    tmpFileStr = strcat(path, file);
    save(tmpFileStr, 'odeSettings', 'odeOut', 'error_runs', 'totalPrcc', ...
        'totalSign', 'totalPrccSignEval', 'paramMatrix', 'icMatrix', ...
        'pulseMatrix', 'LHSmatrix', 'PRCC_var');
end

end

function totalPrccSignEval = prcc_significance_eval(totalPrcc, totalSign, alpha)
%============================================================================    
%
% Evaluate a prcc value against its corresponding significance value and a
% significance threshold value (alpha).
%
% If the signficance is less than the alpha value then the prcc value is
% significant. In this case return the prcc value.
%
% If the signficance is greater than or equal to the alpha value then the
% prcc value is not significant. In this case return NaN (not a number).
%
% If the significance is not a finite number (i.e. it is NaN or Inf) then
% also return NaN.
%
%============================================================================

totalPrccSignEval = totalPrcc;
dims = size(totalPrcc);
for i = 1:dims(1)
    for j = 1:dims(2)
        for k = 1:dims(3)
            
            prcc = totalPrcc(i,j,k);
            sign = totalSign(i,j,k);
            
            if ~isfinite(sign) || sign >= alpha
                totalPrccSignEval(i,j,k) = NaN;
            end

        end
    end
end    
end

function odeOut = initVars(outputCount, timePointCount, runCount)
%============================================================================
%
% Pre-allocate the output cell array for speed.
%
%
% Inputs:
%
%   outputCount: The number of model outputs, the same as the number of model
%                equations, and model initial conditions.
%
%    timePointCount: The number of model output time points for each model run.
%
%   runCount: The number of model runs.
%
% Output:
%
%   odeOut: A cell array that has 1 element per model output. Each element for
%           a model output is a matrix with timePointCount rows and run count
%           columns. Column r is the model results for that model output for
%           all time points for run r. Row t is the model results for that
%           model output for for time point t for all runs.
%
%============================================================================

odeOut = cell(outputCount, 1);
for o = 1:outputCount
    odeOut(o) = {zeros(timePointCount, runCount)};
end

end

function odeOut = storeModelRunOutput(odeOut, run, y)

%============================================================================
%
% Store the model output for a single run, when the solver did not quit early -
% model outputs were produced for all time points.
%
% Inputs:
%
%   odeOut: A cell array that holds model output.
%
%    run: Current run number.
%
%   y: The output vector from the solver for the current run.
%
% Output:
%
%  odeOut: Updated with the model output y.


[timePointCount, outputCount] = size(y);

for i = 1:outputCount
    odeOut{i}(:,run) = y(:,i);
end
end

function odeOut = saveError(odeOut, run, y, tspan)
%============================================================================
%
% When an error occurs, fill odeOut with what model output produced by the
% solver, and pad the missing output with the values from the last time point
% for which a value was produced by the solver.
%
% Inputs:
%
%   odeOut: A cell array that holds model output.
%
%    run: Current run number.
%
%   y: The output vector from the solver for the current run.
%
%   tspan: Time vector of all time points for which output is normally
%          produced.
%
% Output:
%
%  odeOut: Updated with the model output y, and padded with missing model
%          output.
%
%============================================================================


% timePointCount is the number of time points for which the solver produced
% output.
%
% outputCount is the number of model otuput variables (number of model
% equations).

[timePointCount, outputCount] = size(y);
tspanCount = length(tspan);
padCount = tspanCount - timePointCount;
padOnes = ones(padCount);

for i = 1:outputCount
    % Use the value of the output for the last time point a value was produced
    % by the solver.
    pad = padOnes * y(:,end);
    odeOut{i}(:,run) = [y(:,i); pad];
end

end

