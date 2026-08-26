clear; clc; close all;

% This code:
%   1. Computes the equilibrium branch
%   2. Uses a known Hopf location
%   3. Starts a periodic-solution branch from that Hopf
%   4. Continues the periodic branch
%   5. Plots the cycle min/max algae range

%% Add DDE-BIFTOOL

addpath(genpath('/Users/milanarado/Documents/dde_biftool_v3.1.1'));

%% Setup

ind_K = 1;

funcs = set_funcs( ...
    'sys_rhs', @Daphnia_RHS, ...
    'sys_tau', @() 1);

%% Start at known positive equilibrium

K0 = 0.05;
par = K0;

x0 = analytic_steady_guess(K0);

stst.kind = 'stst';
stst.parameter = par;
stst.x = x0;

method = df_mthod(funcs, 'stst');

[stst, success] = p_correc(funcs, stst, [], [], method.point);

if ~success
    error('Initial steady-state correction failed.');
end

%% Build equilibrium branch

branch = df_brnch(funcs, ind_K, 'stst');

branch.parameter.min_bound = [ind_K 0.05];
branch.parameter.max_bound = [ind_K 3.5];

% Larger values are faster
branch.parameter.max_step = [ind_K 0.01];

branch.point = stst;

%% Create second equilibrium point

stst2 = stst;
stst2.parameter(ind_K) = K0 + 0.01;

[stst2, success] = p_correc(funcs, stst2, [], [], method.point);

if ~success
    error('Second steady-state correction failed.');
end

branch.point(2) = stst2;

%% Continue equilibrium branch

branch.method.continuation.plot = 0;

[branch, ~, ~, ~] = br_contn(funcs, branch, 400);

%% Extract equilibrium branch values

K_vals = arrayfun(@(p) p.parameter(ind_K), branch.point);
F_vals = arrayfun(@(p) p.x(1), branch.point);

% Dimensional algae density
F_dim = K_vals .* F_vals;

%% Use known Hopf location

hopf_K_target = 0.78; % guess from other code seeing first stability change

K_branch = arrayfun(@(p) p.parameter(ind_K), branch.point);
[~, hopf_ind] = min(abs(K_branch - hopf_K_target)); % find the K value in the branch that is closest to target, only record K value not distance

fprintf('\nUsing branch point closest to Hopf K = %.10f\n', hopf_K_target);
fprintf('Selected branch index = %d, K = %.10f\n', ...
        hopf_ind, branch.point(hopf_ind).parameter(ind_K));

%% Start periodic-solution branch from Hopf

degree = 3; % cubic polynomial approximation
intervals = 40; % divide a cycle into 20 parts to approximate
radius = 1e-3; % starts the initial cycle very small around the hopf

% Biftool function that creates oscillation around hopf_ind
[psol_branch, suc] = SetupPsol( ...
    funcs, branch, hopf_ind, ...
    'contpar', ind_K, ...
    'degree', degree, ...
    'intervals', intervals, ...
    'radius', radius);

if ~suc
    error('SetupPsol failed.');
end

%% Continue periodic branch forward

psol_branch.parameter.min_bound = [ind_K 0.05];
psol_branch.parameter.max_bound = [ind_K 3.5];

% Larger step is faster
psol_branch.parameter.max_step = [ind_K 0.01];

psol_branch.method.continuation.plot = 0;

[psol_forward, ~, ~, ~] = br_contn(funcs, psol_branch, 400); % predicts the next periodic orbit, then correct it to satisfy DDE = 0 for full period
% just records the new branch value which contains period, k value, shape..
% extrapolates new period, shape, from two previous then newton corrects

%% Continue periodic branch backward

psol_backward = br_rvers(psol_branch);

psol_backward.parameter.min_bound = [ind_K 0.05];
psol_backward.parameter.max_bound = [ind_K 3.5];
psol_backward.parameter.max_step  = [ind_K 0.01];

psol_backward.method.continuation.plot = 0;

[psol_backward, ~, ~, ~] = br_contn(funcs, psol_backward, 200);

%% Combine forward and backward branches

psol_all = psol_forward;
psol_all.point = [psol_backward.point(end:-1:1), psol_forward.point]; % makes one unified branch, helpful for folding when goes backward past hopf

%% Extract periodic-cycle algae ranges

K_psol = nan(size(psol_all.point));
Fmin_psol = nan(size(psol_all.point));
Fmax_psol = nan(size(psol_all.point));

for i = 1:length(psol_all.point)

    p = psol_all.point(i);

    K = p.parameter(ind_K);

    F_profile = p.profile(1,:); % gets algae vals for one period cycle around each point

    F_profile_dim = K * F_profile; % Redimensionalize algae density

    K_psol(i) = K; % store K vals
    Fmin_psol(i) = min(F_profile_dim); % find min of that period
    Fmax_psol(i) = max(F_profile_dim); % find max of that period
end


%% Compute periodic-orbit amplitude

amp_psol = Fmax_psol - Fmin_psol; % difference of max and min gives amplitude for each branch point i

% Split into small- and large-amplitude cycles
amp_split = 0.12;   % arbitrary value to threshold into small and large amplitude cycles

is_small = amp_psol < amp_split;
is_large = amp_psol >= amp_split;

%% Plot Figure 1a-style branches

figure; hold on;

% Equilibrium branch
plot(K_vals, F_dim, 'k:', 'LineWidth', 1.5, ...
    'DisplayName', 'Equilibrium branch');

plot(K_psol(is_small), Fmin_psol(is_small), 'b.', 'MarkerSize', 10, ...
    'DisplayName', 'Small-amplitude cycles');
plot(K_psol(is_small), Fmax_psol(is_small), 'b.', 'MarkerSize', 10, ...
    'HandleVisibility', 'off');

plot(K_psol(is_large), Fmin_psol(is_large), 'r.', 'MarkerSize', 10, ...
    'DisplayName', 'Large-amplitude cycles');
plot(K_psol(is_large), Fmax_psol(is_large), 'r.', 'MarkerSize', 10, ...
    'HandleVisibility', 'off');

xlabel('K (mg C L^{-1})');
ylabel('Algal density (mg C L^{-1})');
title('Algal Density Cycles Around Hopf');
xlim([0.05 3.5]);
set(gca,'YScale','log');   % convert to log scale on y axis
legend('Location','best');
box on;

%% Function


function x_guess = analytic_steady_guess(K)

    q     = 1.0;
    IJ    = 5.23e-3;
    IA    = 1.91e-2;
    fh    = 0.164;
    chi   = 0.77;
    gamma = 1.51e-3;
    sigJ  = 0.49;
    sigA  = 0.43;
    muJ   = 1.09e-5;
    muA   = 5.92e-6;
    w     = 4.8e-3;

    lambda = fh / K;

    alphaJ = (chi / gamma) * (IJ * sigA / q);
    alphaA = (chi / gamma) * (IA * sigA / q);

    etaJ = muJ / (q * IJ * sigJ);
    etaA = muA / (q * IA * sigA);

    rho = (IJ * sigJ) / (w * q);

    zeta = 1000;

    steady_eq = @(F) alphaA*exp(-(etaJ/rho)*((F+lambda)/F)^2) ...
                     - etaA*((F+lambda)/F)^2;

    F_eq = fzero(steady_eq, [1e-5, 0.99]);

    H_eq = ((F_eq + lambda)/F_eq)^2;
    S_eq = exp(-(etaJ/rho)*H_eq);

    B = alphaJ*(1 - S_eq)/(etaJ*H_eq);
    R = (1 - F_eq)*(F_eq + lambda);

    A_eq = R/(1 + B);
    J_eq = B*A_eq;

    Q_eq = H_eq/zeta;

    x_guess = [
        F_eq;
        J_eq;
        A_eq;
        Q_eq;
        Q_eq;
        Q_eq;
        Q_eq;
        Q_eq
    ];
end