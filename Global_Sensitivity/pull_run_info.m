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

function [run_info] = pull_run_info(filename)  
    
    raw = readtable(filename,'sheet','run_info');

    % convert to correct format
    run_info = {};
    for i = 2:size(raw,2)
    
        run_info{end+1} = raw{1,i};
                
    end
        
       
end
    