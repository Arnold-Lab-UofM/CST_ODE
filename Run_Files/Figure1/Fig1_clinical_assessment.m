%% 
fn = '../workspaces/2022-04-29_VALENCIA-16s-to-SS-type.mat';

[multi_stable_pid,multi_stable_types,mono_stable_pid,...
    mono_stable_types] = analyze_Clinical_Transitions(fn);

