%% symbolic_solns.m
% [S, Jmat] = symbolic_solns(num_sp,Type)
%
% Use this to calculate steady states and jacobian gLV
% Input:
%   * N = numbers of species
%
% Output: 
%   * S = calculated steady-states
%   * Jmat = calculated Jacobian
%
% Other Info
% 2) GENERATE EQUATIONS TOT SOLVE FOR SS
% Since the gLV equations exhibit a pattern, we can populate the equations
% within a for loop (allows this to generalize to an "N" species model).
%    
% 3) CALCULATE JACOBIAN TO DETERMINE STABILITY
% Apply linear stability analysis to determine the behavior around the
% equilibrium points. To do this, we need to now determine the Jacobian.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% v2: Feb 19, 2021
% v3: Jan 2, 2023 (changed input variable name)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%
function [S, Jmat] = symbolic_solns(N)

    clear eqns sub

    st_var = sym('y',[1 N]); % define state variables (the species)
    g_var = sym('M',[1 N]); % define growth rates (on per species)
    B_var = transpose(sym('B',[N N])); % define interaction coefficients (species x species matrix)

    for i = 1:length(st_var)
        sub(i) = g_var(i) + B_var(i,i)*st_var(i);
        for j = 1:length(st_var)
            if i ~= j
                sub(i) = sub(i) + B_var(i,j)*st_var(j); % populate interaction term
            end
        end
        eqns(i) = sub(i)*st_var(i) == 0; % set equal to zero to solve for steady states
    end
    S = solve(eqns,st_var); % solves the system of ODEs set to 0
    
    % Jacobian
    clear eqns sub
    for i = 1:length(st_var)
        sub(i) = g_var(i) + B_var(i,i)*st_var(i);
        for j = 1:length(st_var)
            if i ~= j
                sub(i) = sub(i) + B_var(i,j)*st_var(j); % populate interaction term
            end
        end
        eqns(i) = sub(i)*st_var(i);
    end

    Jmat = jacobian(eqns,[st_var]); % analytical jacobian

end