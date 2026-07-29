clear; clc; close all;

%% Analytic-equilibrium + finite-perturbation stability scan
% not linear stability

%% Settings

q     = 1.0;
IJ    = 5.23e-3;
IA    = 1.91e-2;

K_vals = linspace(0.05, 3.5, 120);

perturb = 1e-2; % perturb F by 1 percent, later multiplies current F by 1.01
phi_end = 250;
perturb_var = 2; % 1 is F, 2 is J, 3 is A

state_names = {'F algae', 'J juvenile Daphnia', 'A adult Daphnia', 'Q1', 'Q2', 'Q3', 'Q4', 'Q5'};

opts = ddeset('RelTol', 1e-6, 'AbsTol', 1e-8); % standard accuracy settings

%% Storage

Eq_dim = nan(size(K_vals)); % stores equilibrium values
finite_stable = false(size(K_vals)); % array of falses which will be changed to true if returns to equilibrium
dev_ratio = nan(size(K_vals)); % stores ratios of final/initial deviation
%% Main loop over K

fprintf('\nAnalytic equilibrium finite-perturbation scan\n');
fprintf('Perturbation size = %.1e\n\n', perturb);

for i = 1:length(K_vals)

    K = K_vals(i);

    fprintf('Testing K = %.4f (%d of %d)\n', K, i, length(K_vals));

    %% Compute analytic equilibrium

    xeq = analytic_steady_guess(K); % calculates the equilibrium from analytic equilibrium equation

    % extracts the equilibrium value of chosen variable and dimensionalizes it
    if perturb_var == 1
        Eq_dim(i) = K * xeq(1);          % dimensional algae density
    elseif perturb_var == 2
        Eq_dim(i) = (q*K/IJ) * xeq(2);   % dimensional juvenile Daphnia density
    elseif perturb_var == 3
        Eq_dim(i) = (q*K/IA) * xeq(3);   % dimensional adult Daphnia density
    else
        Eq_dim(i) = xeq(perturb_var);
    end
    %% Apply finite perturbation to algae F

    history = xeq; % set the history of DDE equal to the calculated equilibrium
    history(perturb_var) = xeq(perturb_var) * (1 + perturb); % perturbs F value of history by 1.01

    %% Simulate DDE

    sol = dde23(@(phi,y,Z) rhs_for_dde23(phi,y,Z,K), 1, history, [0 phi_end], opts); % simulate for perturbed history

    phi = sol.x; % time points 
    Y = sol.y; % simulated vals of variables

    %% Distance from equilibrium

    dev = vecnorm(Y - xeq, 2, 1); % calculates the distance between the perturbed state and equilibrium at every phi value using 2-norm Euclidean

    start_dev = median(dev(phi <= 30)); % median to be more resistant to noise, gets initial perturb distance value for below 30
    end_dev   = median(dev(phi >= 150)); % gets final perturb distance for above 150

    dev_ratio(i) = end_dev / start_dev; % calculates the ratio, <1 means got smaller, >1 got bigger

    %% Finite-perturbation stability rule
    % Stable means the perturbation became much smaller by the end, not linear stability tho.

    finite_stable(i) = dev_ratio(i) < 0.1; % stable if perturbation shrank to 10 percent of original
end

%% Plot analytic equilibrium branch with finite-perturbation classification

figure; hold on;

plot(K_vals, Eq_dim, 'w', 'LineWidth', 1.2);

plot(K_vals(finite_stable), Eq_dim(finite_stable), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5); % plots points where we had true for stable array

plot(K_vals(~finite_stable), Eq_dim(~finite_stable), 'ro', 'MarkerSize', 5); % plots points where we had false for stable array

xlabel('K (mg C L^{-1})');
ylabel(sprintf('Equilibrium %s', state_names{perturb_var}));
title(sprintf('%s perturbation %.0e', state_names{perturb_var}, perturb));
xlim([0.05 3.5]);

legend('Analytic equilibrium branch', 'Returns after finite perturbation', 'Does not return after finite perturbation', 'Location', 'best');

%% Plot deviation ratio

figure;
semilogy(K_vals, dev_ratio, 'bo-');
hold on;
yline(1, '--'); % line defining growth or shrink based on which side
yline(0.1, ':'); % line defining return threshold

xlabel('K');
ylabel('Final deviation / initial deviation');
title(sprintf('Finite-perturbation response, perturb = %.0e', perturb));
xlim([0.05 3.5]);

legend('Deviation ratio', 'Line of no response', 'Return threshold', 'Location', 'best');

%% Print rough transition location

unstable_inds = find(~finite_stable); % puts all unstable equilibrium K values

if isempty(unstable_inds)
    fprintf('\nAll tested K values returned after perturbation %.1e.\n', perturb);
else
    fprintf('\nFirst K value that did not return: %.6f\n', K_vals(unstable_inds(1))); % prints the first of those unstable equilibrium values
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dydphi = rhs_for_dde23(phi, y, Z, K)
    xx = [y, Z(:,1)];
    dydphi = Daphnia_RHS_local(xx, K);
end

function dydphi = Daphnia_RHS_local(xx, K)
% Same 8-variable transformed fixed-delay Daphnia-algae model.
%
% xx(:,1) = current state
% xx(:,2) = delayed state
%
% State:
%   F, J, A, Q1, Q2, Q3, Q4, Q5

    y = xx(:,1);
    Z = xx(:,2);

    F  = y(1);
    J  = y(2);
    A  = y(3);
    Q1 = y(4);
    Q2 = y(5);
    Q3 = y(6);
    Q4 = y(7);
    Q5 = y(8);

    A_delay  = Z(3);
    Q1_delay = Z(4);
    Q2_delay = Z(5);
    Q3_delay = Z(6);
    Q4_delay = Z(7);
    Q5_delay = Z(8);

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

    nu   = 1;
    zeta = 1000;

    if F <= 0
        F = eps;
    end

    H = ((F + lambda)/F)^2;

    S = exp(-zeta*(etaJ/rho) * ( ...
          Q1 - Q1_delay ...
        + nu   * (Q2 - Q2_delay) ...
        + nu^2 * (Q3 - Q3_delay) ...
        + nu^3 * (Q4 - Q4_delay) ...
        + nu^4 * (Q5 - Q5_delay) ...
        + nu^5 * Q5));

    dF = (1/rho)*((1 - F)*(F + lambda) - J - A);

    dJ = (alphaJ/rho)*(A - A_delay*S) ...
         - (etaJ/rho)*H*J;

    dA = (alphaA/rho)*A_delay*S ...
         - (etaA/rho)*H*A;

    dQ1 = -nu*Q1 + (1/zeta)*H;
    dQ2 = -nu*Q2 + Q1;
    dQ3 = -nu*Q3 + Q2;
    dQ4 = -nu*Q4 + Q3;
    dQ5 = -nu*Q5 + Q4;

    dydphi = [
        dF;
        dJ;
        dA;
        dQ1;
        dQ2;
        dQ3;
        dQ4;
        dQ5
    ];
end

function x_guess = analytic_steady_guess(K)
% Analytic positive equilibrium of the transformed model.

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