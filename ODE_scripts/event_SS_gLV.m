%% event_SS_probiotic_gLV.m
% Event function for ODE solver to stope the model when it reaches
% steady-state
%
% Will need to modify dy function for each new model
% May need to modify the slope_threshold depending how close to steady
% state you want the model to terminate at 

function [x, isterm, dir] = event_SS_gLV(t,y,params)
% Supply the differential equations
dy = lhs_ode_gLV(t,y,params);

% Determines when the derivative of the solution is close to zero (steady
% state)
slope_threshold = 1E-4; % May vary depending on model (check visually)
x = norm(dy) - slope_threshold;
isterm = 1; % tells function to terminate solver
dir = 0; % directionality
end