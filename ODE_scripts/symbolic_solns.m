%% symbolic_solns.m
% Use this to calculate steady states and jacobian of 'GUASE'-type or 'CAlpha' type
% generalized Lotka Volterra Equations
% [S, Jmat] = symbolic_solns(num_sp,Type)
% Input:
% num_sp = numbers of species
% Type = 'GAUSE' or 'CAlpha' dependence on type of gLV formulation
%
% Output: 
% S = calculated steady-states
% Jmat = calculated Jacobian
%
% Questions can be directed at Christina Lee (chyylee@umich.edu)
% 2/19/2021   
%
% Other Info
% 2) GENERATE EQUATIONS TOT SOLVE FOR SS
% Since the gLV equations exhibit a pattern, we can populate the equations
% within a for loop (allows this to generalize to an "N" species model).
%    
% 3) CALCULATE JACOBIAN TO DETERMINE STABILITY
% Apply linear stability analysis to determine the behavior around the
% equilibrium points. To do this, we need to now determine the Jacobian.

function [S, Jmat] = symbolic_solns(num_sp,Type)
    N = num_sp; % enter number of species in the model, here we want 3


    clear eqns
    clear sub
    if contains(lower(Type),'gause')
        st_var = sym('y',[1 N]); % define state variables (the species)
        g_var = sym('M',[1 N]); % define growth rates (on per species)
        B_var = transpose(sym('B',[N N])); % define interaction coefficients (species x species matrix)
        K_var = sym('K',[1 N]); % define carrying capacities(on per species);
        for i = 1:length(st_var)
            sub(i) = K_var(i) - st_var(i);
            for j = 1:length(st_var)
                if i ~= j
                    sub(i) = sub(i) - B_var(i,j)*st_var(j); % populate interaction term
                end
            end
            eqns(i) = g_var(i)*st_var(i)*(sub(i)/K_var(i)) == 0; % set equal to zero to solve for steady states
        end
        S = solve(eqns,st_var); % solves the system of ODEs set to 0
        
        % Jacobian
        clear eqns sub
        for i = 1:length(st_var)
            sub(i) = K_var(i) - st_var(i);
            for j = 1:length(st_var)
                if i ~= j
                    sub(i) = sub(i) - B_var(i,j)*st_var(j); % populate interaction term
                end
            end
            eqns(i) = g_var(i)*st_var(i)*(sub(i)/K_var(i));
        end

        Jmat = jacobian(eqns,[st_var]); % analytical jacobian
    elseif contains(lower(Type),'calpha')
        st_var = sym('y',[1 N]); % define state variables (the species)
        g_var = sym('M',[1 N]); % define growth rates (on per species)
        B_var = transpose(sym('B',[N N])); % define interaction coefficients (species x species matrix)
        K_var = sym('K',[1 N]); % define carrying capacities(on per species);

        for i = 1:length(st_var)
            sub(i) = g_var(i) - (g_var(i)/K_var(i))*st_var(i);
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
            sub(i) = g_var(i) - (g_var(i)/K_var(i))*st_var(i);
            for j = 1:length(st_var)
                if i ~= j
                    sub(i) = sub(i) + B_var(i,j)*st_var(j); % populate interaction term
                end
            end
            eqns(i) = sub(i)*st_var(i);
        end

        Jmat = jacobian(eqns,[st_var]); % analytical jacobian
        
    elseif contains(lower(Type),'ralpha')
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
    else
        disp(['Please type: GAUSE, calpha, Ralpha'])
    end

end