%% G1_Generate_Reference_Virtual_Patient.m
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Goal: Use the initial/base sensitivity analysis to generate a larger
% reference population that can be used to match population-level
% observations in equilibrium behaviors (EB).
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Description of code workflow:
%   (1) Load the base/initial workspace
%           - Input: 'SSConfig-Analysis-XXXX.mat'
%   (2) Use parameter sets of known EB to generate new parameter sets
%           - Input: LHSmat selected by EB type
%           - Output: Larger parameter matrix with ("Reference"
%           population")
%   (3) Match the clinically observed EB Frequencies
%           - Input: EB frequencies for HMP or Gajer Cohort (use result of
%               G1_Analyze_Clinical_Equilibrium_Behaviors.m)
%           - Output: New SSConfig-Analysis.mat that has the matched
%           population data [WARNING: may overwrite previous workpaces]
%
% REQUIRES: [Add folders to MATLAB path]
%       * analyze_global_SS.m
%       * calc_SS_stability.m
%           - requires "S" which is the models symbolic equations, "Jmat" which
%           is the Jacobian of the system of ODEs
%       * analyze_Global_CST_SS.m
%   * MATLAB TOOLBOXES: Parallel Computing Toolbox, Bioinformatics Toolbox,
%           Symbolic Math Toolbox
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee (chyylee@umich.edu)
% University of Michigan
% Apr 4, 2023 (Created during revision process)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% 1. Loads required workspace and sets up base variables
clear;clc;
fdr_loc = '../workspaces/'; 
load(strcat(fdr_loc,'SSConfig-Analysis-Model_LHS_10x.mat'))

if ~exist('SS_names_CST','var') % checks for variable name change across versions
    SS_names_CST = SS_namesv;
    all_nm_CST = all_nm;
end

EB_names = {'1SS: [Li] CST-III';
    '1SS: [oLB] CST-I/II/V';
    '1SS: [NO] CST-IV';
    '2SS: [NO] CST-IV or [Li] CST-III';
    '2SS: [NO] CST-IV or [Li] CST-III';
    '2SS: [NO] CST-IV or [Li] CST-III'};

%% 2. Iterate through each EB and create 5000 samples
NR = 5000; % Number of "samples", our work uses 5000
rng(1);
full_paramMatrix = NaN(NR*length(EB_names),length(param_names));
c = 1;
for indx = 1:length(length(EB_names))
    [sel_idx] = find(all_nm_CST == EB_names(indx)); % all_nm has been updated to all_nm_CST in some versions of code (Number of LHS samples x 1 string array)
    sel_nets = LHSmat(sel_idx,:); %final parameter sets to use
   
    % 2.) Prep for LHS
    LHSmin = mean(sel_nets);
    LHSmax = std(sel_nets);
    
    % format for LHS sampling
    params = {};
    for i = 1:length(LHSmin)
        tmp = {param_names{i}, 'n',LHSmin(i),LHSmax(i)};
        params{end+1} = tmp;
    end

    paramMatrix = defineMatrix(params, NR,'parameter'); % Generates LHS parameter sets
    full_paramMatrix(c:c+NR-1,:) = paramMatrix;
    c = c + NR;
end

%% 3. MATCH THE POPULATION TO CLINICAL FREQUENCIES
ub = 0.99; % Upper bound of the dominant species
lb = 0.6; % Lower bound of the dominant species
NumPats = 1200; % Total size of the simulated populatiion

% ~~~~~~~~~~ Gajer EB Frequencies ~~~~~~~~~~ (use result of
% G1_Analyze_Clinical_Equilibrium_Behaviors.m)
gajer32 = [0.25 0.3125 0.2188 0.1562,0,0.0625];
counts_per_EB = ceil(NumPats.*gajer32); % 1SS nAB, 1SS Li, 1SS oLB, 2SS Li or nAB, 2SS Li or oLB, 2SS oLB or nAB
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ~~~~~~~~~~ HMP EB Frequencies ~~~~~~~~~~ (uncomment to run HMP)
% hmp101 = [0.307 0.109 0.337 0.158,0.04,0.05];
% counts_per_EB = ceil(NumPats.*hmp101); % 1SS nAB, 1SS Li, 1SS oLB, 2SS Li or nAB, 2SS Li or oLB, 2SS oLB or nAB
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sel_cst_names = {'1SS: [NO] CST-IV';
    '1SS: [Li] CST-III';
    '1SS: [oLB] CST-I/II/V';
    '2SS: [NO] CST-IV or [Li] CST-III';
    '2SS: [Li] CST-III or [oLB] CST-I/II/V';
    '2SS: [NO] CST-IV or [oLB] CST-I/II/V'};

% Remove set with negative growth rates:
neg_gr = sum(LHSmat(:,1:3) < 0,2) > 0;
pos_si = sum(LHSmat(:,[4 8 12]) > 0,2) > 0;

s = RandStream('mlfg6331_64'); 
c = 1;

% describes indexes of species that are dominant (more than one 1 in a row
% indicates multi-stabiilty)
sp_combos = [1 0 0;
    0 1 0;
    0 0 1;
    1 1 0;
    0 1 1;
    1 0 1];

final_population = NaN(NumPats,length(param_names));
for cst_id = 1:length(sel_cst_names)
    indx = find(all_nm_CST == sel_cst_names{cst_id} & ~neg_gr & ~pos_si);
    tmp = StbleSS(indx);
    abundance = cell2mat([tmp{:}]');
    rel_abundance = abundance ./ sum(abundance,2);

    % check which species can be dominant, check lower and upper bounds
    ssu = find(sp_combos(cst_id,:));
    spComp = NaN(size(tmp,1),length(ssu));
    ssn = repmat(ssu,size(tmp,1),1);
    for i = 1:length(ssu)
        rtmp = rel_abundance(ssn==ssu(i),ssu(i));
        bid = rtmp > lb & rtmp < ub;
        spComp(:,i) = bid;
    end
    poss_indx = sum(spComp,2) == 1;

    % Sample from the sets that meet the lower and upper bounds
    num_samples = counts_per_EB(cst_id);
    tmp_id = randsample(s,indx(poss_indx),num_samples);
    final_population(c:c+num_samples-1,:) = LHSmat(tmp_id,:);
    c = c + num_samples;
end

analyze_Global_SS(final_population,{'nAB','Li','oLB'}) % Call function that determines SS configurations


