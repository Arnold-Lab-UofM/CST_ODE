%% lhs_ode_settings_GUI.m

% Modification to Kirschner Group's code to allow for GUI directed input
% for LHS-PRCC analysis. This code requires input using a spreadsheet that
% has information of parameter and initial conditions probability
% distributions.

% Enter spreadsheet in format:
% Sheet 1 (labeled "parameters")
%   First row: column headers (parameter name, distribution, value 1, value 2)
%   Second - end rows: parameters

% Sheet 2 (labeled "initial_conditions")
%   First row: column headers (name, distribution, value 1, value 2)
%   Second - end rows: initial conditions for each ode 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% Jan 12, 2020
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function  ode = lhs_ode_settings_ng()
    dims = [1 50];
    prompt = {'Name of spreadsheet:'};
    dlgtitle = 'Enter LHS datasheet info';
    definput = {'lhs_settings_input.xlsx'};
    flnm = char(inputdlg(prompt,dlgtitle,dims,definput));
    
    [run_info] = pull_run_info(filename);
    
% This is an example of a new format ODE LHS settings file.

% To use this file for a new model do the following:
%     Copy this file and edit it to define the items as needed for the new
%     model.  That includes deleting these header comments and replacing them
%     with ones that describe the new model.
%
%     Create a Matlab file to contain the equations,  Be sure to define item
%     ode.odeModelHandle in this file accordingly.



% To run this model, use the following command from within a Matlab command
% window:
% lhs_ode_run_new('name_of_this_file')

% You can edit the the number of runs, time points to analyze, add or change
% ranges for the parameters, initial conditions and pulse items, as described
% in the comments below. 

% In this simple model all parameters are varied according to a normal
% distribution ('n' in the parameter items). Our production models will
% typically not vary all parameters and most tend to be varied according to a
% uniform distribution ('u').

% This settings file does not vary the initial conditions. You can experiment
% with varying initial conditions by specifying a distribution to use (ex. 'u'
% instead of '') and values that specify a range (uniform 'u' or log uniform
% 'lu' distribution) or that specify a mean and standard deviation (normal
% distribution, 'n').


%===============================================================================
%
% Output:
%   The function return value is a structure named ode, with the following
%   members, which need to be define in this function.
%
%   NR
%       The number of times to run the ode model.
%
%   tspan
%       A vector of time points to produce model run results for.
%
%       This is an input to the ODE solver.
%ode.analysisTimePoints=
%   analysisTimePoints
%       A vector of time points to perform an analysis for.
%
%       Each element of this vector must also be an element of tspan.
%
%   odeModelHandle
%       The name of the Matlab function that contains the model equations.
%       This is an input to the ODE solver.
%       This is required.
%          
%   computeParamsInitCondHandle
%       The name of a Matlab function for computing some parameter values
%       and/or initial conditions based on other parameter values and/or
%       initial conditions.
%
%       This is not required. Leave undefined if not needed.
%
%  solverTimeLimit
%      The amount of time to allot for solving one set of parameters and initial
%      conditions. If the Matlab ODE solver takes longer than this amount the
%      solve is stopped and noted in an error variable.
%
%      In seconds. Can be fractional, ex. 1.5 or 0.5.
%
%  saveOnError
%      Whether or not to save the ODE solution results when an error occurs
%      during a call to the Matlab ODE solver.
%      1 = yes
%      0 = no.
%
%      If yes, pad the output matrix for the solution time steps that were not
%      defined, so the result data structures for all parameter sets will have
%      the same size and shape.Model_LHS.mat
%
%      If no, a failed solver call results in that run's results not being
%      saved in the output matrix.
%
%   parameters
%       A cell array for specifying the parameters.
%
%       There is one element of the cell array for each parameter.
%       The element for a parameter is also a cell array that has the
%       parameter's name, probability distribution and range.
%
%       The probability distribution is one of the following:
%           ''    No distribution, use the first value in the range.
%                 Both values in the range must be the same.
%
%           'u'    Use the 2 values in the range as the min and max for a
%                  uniform distribution.
%                  min must be less than max.
%
%           'lu'   Use the 2 values in the range as the min and max for a
%                  log uniform distribution.
%                  min must be less than max.
%
%           'n'    Use the 2 values in the range as the mean and std.dev.
%                  for a normal distribution.
%                   Std. dev. must not be 0.
%
%       The range always has 2 values, as described above in the probability
%       distribution.
%
%       If a parameter is a computed parameter (will be defined in the function
%       specified by compute_params_handle), then specify a probability
%       distribution of '' and a range of 0.0, 0.0. It doesn't really
%       matter, since any random value defined for it will be overwritten by
%       the computed value, but it would be confusing to specify a range if
%       that range won't be used.
%
%   initialConditions
%       A cell array for specifying the initial conditions.
%
%       This has the same format as the parameter cell array.
%
%       Initial conditions can be chosen from a range, using a probability
%       distribution, or be computed, just as for parameters.
%
%   outputLabels
%       A cell array with string labels for each model output.
%       Usually use the default of 'O1', 'O2', ...
%


% Sample size.
ode.NR = run_info{1,1};

% Time span of the simulation.lhs_ode_settings_predator_prey_new
t_end = run_info{1,2}; % length of the simulations (in days)
ode.tspan = (0:1.0:t_end); % time points where the output is calculated
ode.analysisTimePoints = (0:1.0:t_end); % time points of interest for the uncertainty/sensitivity analysis

% The alpha value to use for PRCC.
% Undefine this (comment out) if no PRCC is to be performed.

% Doing PRCC is generally used for production runs.
%
% Not doing PRCC can be useful for testing and debugging the code which reads a
% settings file like this one and runs the model.
ode.alpha = run_info{1,3};

% A function handle for the Matlab ODE model to analyze.
% For example: odeModelHandle = @ODEmodel;
%
% Be sure to spell this correctly, and that the file containing this
% function exists and is on the Matlab path. 
%
% Don't put it in quotes, ex. odeModelHandle = '@ODEmodel';
% is not correct.
%
ode.odeModelHandle = str2func(strcat('@',run_info{1,4}{1}));

% Optional function handle for computing parameters and/or initial conditions.
% Leave undefined if not needed.
%ode.computeParamsInitCondHandle = @lhs_ode_compute_params_predator_prey

% Redefine the ODE solver time limit, if you want a value different than the                                  
% default. In seconds. Can be fractional, ex. 1.5 or 0.5. 
ode.solverTimeLimit = '60.0';

% Set saveOnError to define how to handle the case where the ODE solver
% fails. saveOnError = 1 means to pad the output matrix for the solution
% time steps that were not defined. saveOnError = 0 means dont' pad the
% oututp matrix.
ode.saveOnError = 0;

% Define the model parameters: {name, distribution, value1, value2}
%
% NOTE:
% Do not vary here any parameters that are based on other parameters and/or
% initial conditions.  Any variation done here will be overwritten when they
% are defined based on those parameters and/or initial conditions.
%
% No blank lines are allowed in between the cell arrays for the individual
% parameters.
%
% Comments lines must be preceded by "...".


[params,ICs] = pull_LHS_parameters(flnm);
ode.parameters = params;
ode.initialConditions = ICs;


% % Labels to use for model outputs.
%
% There should be one for each model equation.
%
% The default labels are 'O1', 'O2', ..., 'Oneq' where neq is the number of
% model equations, which should be the same as the number of initial
% conditions.
icCount = length(ode.initialConditions);
ode.outputLabels = lhs_ode_default_output_labels_new(icCount);


end % RUN: lhs_ode_run_new('filename')