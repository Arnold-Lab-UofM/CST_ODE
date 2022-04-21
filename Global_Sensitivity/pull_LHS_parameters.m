%% pull_LHS_parameters.m
% Function so that the lhs_ode_settings_from_xlsxfile.m can pull parameters
% from an excel spreadsheet. Excel spreadsheets will be useful once we
% start to get a lot of parameters in the system.

% The spreadsheet MUST be in the following format:
% sheets are labeled "parameters" and "initial_conditions".
% the first row of each spreadsheet is a header.
% the columns are in the order of "name of parameters", "distribution", 
% "value 1", "value 2"

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% Jan 8, 2020
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [parameters,initialConditions] = pull_LHS_parameters(filename)  
    
    raw = readtable(filename,'sheet','parameters');

    % convert to correct format
    parameters = {};
    for i = 1:size(raw,1) % first row is header, start at row 2
        nrow = {raw{i,1}{1},raw{i,2}{1},raw{i,3},raw{i,4}};
        if isnan(nrow{2})
            nrow{2} = '';
        end
        parameters{end+1} = nrow;
    end
%%
    raw = readtable(filename,'sheet','initial_conditions');
    initialConditions = {};
    for i = 1:size(raw,1) % first row is header
        if isnan(raw{i,2})
            nrow = {raw{i,1}{1},'',raw{i,3},raw{i,4}};
        else
            nrow = {raw{i,1}{1},raw{i,2}{1},raw{i,3},raw{i,4}};
        end
        initialConditions{end+1} = nrow;
    end
end
    