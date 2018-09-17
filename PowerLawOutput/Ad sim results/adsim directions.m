% Run adsorption simulations

% Concentration input does not change
C = [1e-12,	5e-12,	1e-11,	1e-10,	1e-09,	2e-09,	3e-09,	5e-09,	1e-08]; % M

% MW
MW1 = 66000; % Daltons (BSA)
% MW2 = 430000; % Daltons (FB)

% Model is either 'exp1' or 'gp'
% For 'exp1', choose from 0.1, 1, 10, or 100 as input for params
% This is l = 0.1, 1, 10, and 100, respectively
model = 'exp1';
params = 0.1;

% For 'gp', choose from [-4, 1, 0], [-2, 1, 0], [1, 1, 0], [3, 1, 0]
% This is a = 0.75, 0.5, 2, and 1.2, respectively
% model = 'gp';
% params = [-4, 1, 0]; % a = 0.75

% Sticking probability - glass: 5E-4. Pluronic: 1E-4. Silane: 1E-2.
S = 5E-4;

sigma_l01 = runadsim(C,MW,model,params,S); % Coverage fraction for lambda = 1
