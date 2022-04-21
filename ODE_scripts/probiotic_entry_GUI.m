%% probiotic_entry_GUI.m

% Easier propogation of LHS parameter information.

% GUI will ask for inputs for number of species in the system and how many
% species total are in the system. The last LB species is automatically
% considered the "Probiotic". Likewise could also be considered a "Pathogen
%" in future simulations.

% This version assumes that the parameters can be populated as:
%   1) All BV species are sampled from same probability distribution
%   2) All LB species are sampled from same probability distribution
%   3) Interaction term ranges assume same probability distributions for
%       these classes of interactions (BV on BV, BV on LB, BV on P, LB on BV,
%       LB on LB, LB on P, P on BV, P on LB, P on P)

% Future updates will allow for modifications of subsets of these groups,
% but for now, if you want to edit subpopulations of BV or LB species, it
% can be completed within the excel file manually.

% Output: The name of the excel spreadsheet generated (flnm) and the type
% of species in the system as a vector (sp_type).

% INSTRUCTIONS: Generate input needed for the lhs_ode_settings_GUI or
% lhs_ode_settings_fromxlsx file.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% Jan 12, 2020
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


function [flnm,sp_type] = probiotic_entry_GUI()

%% Ask number of species, how many BV?

    prompt = {'Enter Number of Species:','Enter Number that are BV-associated:'};
    dlgtitle = 'Input';
    dims = [1 35];
    definput = {'4','2'};
    sp_info = inputdlg(prompt,dlgtitle,dims,definput);

    %%
    sp_num = str2double(sp_info{1});
    num_BV = str2double(sp_info{2});
    num_LB = sp_num - num_BV;

    sp_type = cellstr([string(repmat('BV',num_BV,1)); string(repmat('LB',num_LB,1))]);

    % number names
    sp_names = sp_type;
    sp_names{end} = 'P';
    b1 = 1;
    b2 = 1;
    for i = 1:sp_num
        if i <= num_BV
            sp_names{i} = strcat(sp_names{i},num2str(b1));
            b1 = b1 + 1;
        elseif i < sp_num
            sp_names{i} = strcat(sp_names{i},num2str(b2));
            b2 = b2 + 1;
        else
            sp_names{i} = sp_names{i};
        end

    end

    pnames = {'k_{grow}','K'};
    nm_out = {};
    for i = 1:length(pnames)
        for j = 1:length(sp_names)
             temp = strcat(pnames(i),'-',sp_names{j});
             nm_out{end+1} = temp;
        end
    end


    %% Keep growth parameters same across all species types?
    % - All same for BV & all same for LB
    prompt = {'BV Distribution:','BV Value 1:','BV Value 2:',...
        'LB Distribution:','LB Value 1:','LB Value 2:', ...
        'Probiotic Distribution:','Probiotic Value 1:','Probiotic Value 2:'};
    dlgtitle = 'Input Growth Rate';
    dims = [1 50];
    definput = {'u','0.2','1','u','0.2','1','u','0.2','1'};
    gr_info = inputdlg(prompt,dlgtitle,dims,definput);

    BVall = {gr_info{1},str2double(gr_info{2}),str2double(gr_info{3})};
    LBall = {gr_info{4},str2double(gr_info{5}),str2double(gr_info{6})};
    Pinfo = {gr_info{7},str2double(gr_info{8}),str2double(gr_info{9})};

    grow = vertcat(repmat(BVall,num_BV,1),repmat(LBall,num_LB-1,1),...
        Pinfo);

    % Repeat entry process
    prompt = {'BV Distribution:','BV Value 1:','BV Value 2:',...
        'LB Distribution:','LB Value 1:','LB Value 2:',...
        'Probiotic Distribution:','Probiotic Value 1:','Probiotic Value 2:'};
    dlgtitle = 'Input Carrying Capacity';
    dims = [1 50];
    definput = {'lu','1','10','lu','1','10','lu','1','10'};
    gr_info = inputdlg(prompt,dlgtitle,dims,definput);

    BVall = {gr_info{1},str2double(gr_info{2}),str2double(gr_info{3})};
    LBall = {gr_info{4},str2double(gr_info{5}),str2double(gr_info{6})};
    Pinfo = {gr_info{7},str2double(gr_info{8}),str2double(gr_info{9})};

    CCap = vertcat(repmat(BVall,num_BV,1),repmat(LBall,num_LB-1,1),...
        Pinfo);


    %% Fill in matrix;
    dims = [1 100];
    prompt = {'Distribution', 'Value 1', 'Value 2'};
    dlgtitle = 'Input Interaction terms BV on BV';
    definput = {'lu','1E-1','2'};
    BVBV_info = inputdlg(prompt,dlgtitle,dims,definput);

    dlgtitle = 'Input Interaction terms BV on LB';
    definput = {'lu','1E-2','2'};
    BVLB_info = inputdlg(prompt,dlgtitle,dims,definput);
    
    dlgtitle = 'Input Interaction terms BV on P';
    definput = {'lu','1E-2','2'};
    BVP_info = inputdlg(prompt,dlgtitle,dims,definput);

    dlgtitle = 'Input Interaction terms LB on LB';
    definput = {'lu','1E-2','2'};
    LBLB_info = inputdlg(prompt,dlgtitle,dims,definput);

    dlgtitle = 'Input Interaction terms LB on BV';
    definput = {'lu','1E-3','2'};
    LBBV_info = inputdlg(prompt,dlgtitle,dims,definput);
    
    dlgtitle = 'Input Interaction terms LB on P';
    definput = {'lu','1E-3','2'};
    LBP_info = inputdlg(prompt,dlgtitle,dims,definput);
    
    dlgtitle = 'Input Interaction terms P on BV';
    definput = {'lu','1E-3','2'};
    PBV_info = inputdlg(prompt,dlgtitle,dims,definput);
    
    dlgtitle = 'Input Interaction terms P on LB';
    definput = {'lu','1E-3','2'};
    PLB_info = inputdlg(prompt,dlgtitle,dims,definput);
    
    dlgtitle = 'Input Interaction terms P on P';
    definput = {'','1','1'};
    PP_info = inputdlg(prompt,dlgtitle,dims,definput);

    BVBVall = {BVBV_info{1},str2double(BVBV_info{2}),str2double(BVBV_info{3})};
    BVLBall = {BVLB_info{1},str2double(BVLB_info{2}),str2double(BVLB_info{3})};
    BVPall = {BVP_info{1},str2double(BVP_info{2}),str2double(BVP_info{3})};
    LBLBall = {LBLB_info{1},str2double(LBLB_info{2}),str2double(LBLB_info{3})};
    LBBVall = {LBBV_info{1},str2double(LBBV_info{2}),str2double(LBBV_info{3})};
    LBPall = {LBP_info{1},str2double(LBP_info{2}),str2double(LBP_info{3})};
    PBVall = {PBV_info{1},str2double(PBV_info{2}),str2double(PBV_info{3})};
    PLBall = {PLB_info{1},str2double(PLB_info{2}),str2double(PLB_info{3})};
    PPall = {PP_info{1},str2double(PP_info{2}),str2double(PP_info{3})};

    p1 = vertcat(repmat(BVBVall,num_BV,1),repmat(BVLBall,num_LB-1,1),...
        BVPall);
    p2 = vertcat(repmat(LBBVall,num_BV,1),repmat(LBLBall,num_LB-1,1),...
        LBPall);
    p3 = vertcat(repmat(PBVall,num_BV,1), repmat(PLBall,num_LB-1,1),...
        PPall);

    full_int = vertcat(repmat(p1,num_BV,1),repmat(p2,num_LB-1,1),p3);

    %% Initial Conditions

    prompt = {'Distribution', 'Value 1', 'Value 2'};
    dlgtitle = 'Initial Species Abundance';
    definput = {'','2','2'};
    ICsp_info = inputdlg(prompt,dlgtitle,dims,definput);

    spsIC = {ICsp_info{1},str2double(ICsp_info{2}),str2double(ICsp_info{3})};

    full_ICs = vertcat(repmat(spsIC,sp_num,1));

    
    IC_nms = sp_names;

    %% reassemble
    int_names = generate_coeff_labels('fx',sp_names)';

    full_nms = vertcat(nm_out',int_names);
    full_params = vertcat(grow,CCap,full_int);

    hdr = {'name','distribution','value 1','value 2'};

    merge_array = horzcat(full_nms,full_params);
    wtitle_array = vertcat(hdr,merge_array);

    IC_array = horzcat(IC_nms,full_ICs);
    wtitle_ICarray = vertcat(hdr,IC_array);


    %% write to excel sheet

    prompt = {'Name of output file:'};
    dlgtitle = 'Name File';
    definput = {'lhs_settings_input.xlsx'};
    flnm = char(inputdlg(prompt,dlgtitle,dims,definput));

    %%
    writetable(cell2table(merge_array,'VariableNames',hdr),flnm,'sheet','parameters')
    writetable(cell2table(IC_array,'VariableNames',hdr),flnm,'sheet','initial_conditions')


end