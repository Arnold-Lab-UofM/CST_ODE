%% Calculate PRCCs and their significances
%%
%% Inputs:
%%
%% N is the number of runs
%%
%% k is the number of varied parameters, initial conditions.
%%
%% LHSmatrix: LHS matrix (N x k)
%%            The k varied parameter and varied initial condition values for
%%            each of N runs.
%%
%%            LHSmatrix(:,i) are the values for varied parameter/initial
%%            condition i for the N runs.
%%
%% Y: output matrix (time x N), where time is the number of time points the
%%    model was solved for.
%%
%% s: time points to test (row vector)
%%    This is indices into a time span.
%%    If s(i) = l then that means use the lth time point.
%%
%% Outputs:
%%
%% prcc: The prcc values, length(s) x k (time points tested x # of varied
%%       items).
%%
%% prcc_significance: The significance (p-values) of the prcc values,
%%                    length(s) x k (time points tested x # of varied items).
%%

%%
%% Originally by Simeone Marino, May 29 2007
%% N.B.: this version uses ONE output at a time

function [prcc, prcc_significance] = lhs_ode_prcc_new(LHSmatrix, Y, s)

timePointCount = length(s);

% Choose the time points (rows of Y) of interest, and transpose so the rows are
% the runs and the columns are the time points of interest, i.e. change Y(time
% x N) to Y(N, timePointCount).
Y = Y(s,:)';

% k is the number of varied parameters/initial conditions.
% timePtCount is the number of model outputs.
[~, k] = size(LHSmatrix);
[~, timePtCount] = size(Y);
if timePtCount ~= timePointCount
    msg = ['After selecting time points from Y, the number of columns of', ...
           ' Y, %d, is not the same as the length of s, %d'];
    fprintf(msg, timePtCount, timePointCount);
    return;
end

% Calculate PRCCs and significances for each of the k varied parameters/initial
% conditions.
prcc = zeros(k, timePtCount);
prcc_significance = zeros(k, timePtCount);
for i = 1:k
    % The LHS matrix with the column for parameter i deleted.
    z = LHSmatrix;
    z(:,i) = [];

    % LHSmatrix(:,i) is Nx1, z is Nx(k-1) and Y is N x timePointCount
    % rho and p are 1 x timePointCount.
    [rho, p] = partialcorr(LHSmatrix(:,i), Y, z, 'type', 'Spearman');
    
    prcc(i,:) = rho;
    prcc_significance(i, :) = p;
end

% Change prcc and prcc_significance to be timePtCount x k.
prcc = prcc';
prcc_significance = prcc_significance';

