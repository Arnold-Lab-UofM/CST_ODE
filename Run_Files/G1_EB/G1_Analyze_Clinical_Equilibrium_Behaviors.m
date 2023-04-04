%% G1_Analyze_Clinical_Equilibrium_Behaviors
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Goal: Use the initial/base sensitivity analysis to generate a larger
% reference population that can be used to match population-level
% observations in equilibrium behaviors (EB).
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Description of code workflow:
%   (1) Define centroids and names of centroids
%   (2) Run neareast centroid classifier for HMP data
%       - Requires the "2022-04-29-VALENCIA-16s-to-SS-type.mat"
%           - linearTransMap: # of patients x 9 matrix with the
%               frequenciess of each state transition
%                   Columns; (1) CST -IV to CST -IV (2) CST -IV to CST
%                   -III (3) CST -IV to CST -I/II/V (4) CST -III to CST -IV
%                   (5) CST -III to CST -III (6) CST -III to CST -I/II/V
%                   (7) CST -I/II/V to CST -IV (8) CST -I/II/V to CST -III 
%                   (9) CST -I/II/V to CST -I/II/V
%       - Output:
%           - allAssigned: # of patients x 1 array of numerical assignment
%           of EB
%           - summaryData: Table with summary of results
%   (3) Repeat process with Gajer
%
% REQUIRES: get_clinical_EB.m
%   
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee (chyylee@umich.edu)
% University of Michigan
% Apr 4, 2023 (Code documented during revision process)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 1. Input Centroid Values and Associated EB Names

% SSref: Centroids
SSref = [1 0 0;
    0 1 0;
    0 0 1;
    0.5 0.5 0;
    0 0.5 0.5;
    0.5 0 0.5];

CentroidNames = {'1SS: nAB Dominated','1SS: Li Dominated', '1SS: oLB Dominated',...
    '2SS: nAB or Li dominated', '2SS: Li or oLB dominated', ...
    '2SS: nAB or oLB dominated'};

%% 2. HMP Transitions
load('../workspaces/HMPData.mat','linearTransMap',...
 'PatSamNum');

L = 10;
less_10 = PatSamNum < L; %Remove Samples with less than 10 time points
[AssignedState_HMP,summaryData_HHMP] = get_Clinical_EB(linearTransMap(~less_10,:),...
    SSref,CentroidNames);

 %% Gajer Transtions
load('../workspaces/GajerData.mat','linearTransMap');
[AssignedState_Gajer,summaryData_Gajer] = get_Clinical_EB(linearTransMap,...
    SSref,CentroidNames);


