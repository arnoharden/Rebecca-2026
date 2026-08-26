clear; clc; close all;

%% Add DDE-BIFTOOL

addpath(genpath('/Users/milanarado/Documents/dde_biftool_v3.1.1'));

%% Setup

ind_K = 1;

funcs = set_funcs('sys_rhs', @Daphnia_RHS, 'sys_tau', @() 1); % BIFTOOL function that sets other script with fixed delay model as RHS

%% Start at a known positive equilibrium

K0 = 0.05;
par = K0;

x0 = analytic_steady_guess(K0);

stst.kind = 'stst';
stst.parameter = par;
stst.x = x0;

method = df_mthod(funcs, 'stst'); % BIFTOOL function defining method to solve specifically steady states

[stst, success] = p_correc(funcs, stst, [], [], method.point); % BIFTOOL function that adjusts the stst.x so that RHS is zero

if ~success
    error('Initial steady-state correction failed.'); % stops if failed correction
end

% disp(stst.x)

%% Build branch

branch = df_brnch(funcs, ind_K, 'stst'); % creates branch which is a fancy array including equilibrium point for all and parameter value for each (i) for equilibria

branch.parameter.min_bound = [ind_K 0.05];
branch.parameter.max_bound = [ind_K 3.5];
branch.parameter.max_step = [ind_K 0.005];

branch.point = stst;

%% Create second point

stst2 = stst;
stst2.parameter(ind_K) = K0 + 0.005;

[stst2, success] = p_correc(funcs, stst2, [], [], method.point);

if ~success
    error('Second steady-state correction failed.');
end

branch.point(2) = stst2;

%% Disable BIFTOOL's live continuation plot

branch.method.continuation.plot = 0;

%% Continue forward

[branch, ~, ~, ~] = br_contn(funcs, branch, 800); % takes the previous two branch.point entries and guesses the next K value and eq variable solution then corrects to have RHS=0 (Newton)
%% Improve stability computation settings

branch.method.stability.minimal_real_part = -1;
branch.method.stability.max_number_of_eigenvalues = 1000;
branch.method.stability.max_newton_iterations = 100;
branch.method.stability.root_accuracy = 1e-11;

%% Compute stability

branch = br_stabl(funcs, branch, 0, 1); % linearizes near each steady state as function of perturbation then plugs in perturbation and solves for characteristic roots

%% Extract values

K_vals = arrayfun(@(p) p.parameter(ind_K), branch.point);
F_vals = arrayfun(@(p) p.x(1), branch.point);
J_vals = arrayfun(@(p) p.x(2), branch.point);
A_vals = arrayfun(@(p) p.x(3), branch.point);

F_dim = K_vals .* F_vals; % re-dimensionalize by reversing the original non-dim equation

%% Count unstable characteristic roots
% chracteristic roots are small perturbations around the equilibrium
% point

nunst = zeros(size(branch.point));
max_real = zeros(size(branch.point));

for i = 1:length(branch.point)
    roots_i = branch.point(i).stability.l0; % where br_stabl stored the characteristic roots

    real_parts = real(roots_i);

    nunst(i) = sum(real_parts > 1e-10); % count the number of unstable roots bc realpart > 0 means unstable, 0.000001 to avoid little noise around 0
    max_real(i) = max(real_parts); % creates array of the max real part of each set of characteristic roots for each branch point
end

%% Local dde23 wrapper for your existing Daphnia_RHS.m

function dydphi = rhs_for_dde23(phi, y, Z, par)
xx = [y, Z(:,1)];
dydphi = Daphnia_RHS(xx, par);
end
%% Plot steady-state branch in Figure 1a variables

figure;
plot(K_vals, F_dim, 'b.', 'MarkerSize', 9);
xlabel('K (mg C L^{-1})');
ylabel('Equilibrium algal density (mg C L^{-1})');
title('BIFTOOL steady-state branch');
xlim([0.05 3.5]);

%% Plot max real characteristic root
% Shows when steady state is stable or unstable since it is only unstable
% when real of all characteristic roots < 0, so if max real > 0 we know
% unstable and if max real < 0 we know stable.

figure;
plot(K_vals, max_real, 'b.-');
hold on;
yline(0, '--');
xlabel('K');
ylabel('Max real part of characteristic roots');
title('Steady-state stability from BIFTOOL');
xlim([0.05 3.5]);

%% Plot number of unstable roots
% Often complex conjugate pairs which is why we see even values of unstable
% roots for each branch point.

figure;
plot(K_vals, nunst, 'o-');
xlabel('K');
ylabel('Number of unstable characteristic roots');
title('Detected stability changes');
xlim([0.05 3.5]);


%% Helper function
% This is just used for the very first equilibrium guess, afterwords the
% continuations use the previous branch points corrected equilibrium as the
% guess.

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

    steady_eq = @(F) alphaA*exp(-(etaJ/rho)*((F+lambda)/F)^2) - etaA*((F+lambda)/F)^2;

    F_eq = fzero(steady_eq, [1e-5, 0.99]); % find value of F in range where steady_eq(F) = 0

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
