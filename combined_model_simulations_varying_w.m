%% Creating Tomas-esque Figures Using Combined Stoichiometric Model
% Constant K with three developmental-threshold values
%
% Daphnia = black
% Algae = dashed black
% Adult proportion = red
% Cycle period = dashed gray
% Stage duration = solid black

clear;
close all;
clc;

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

%% Parameter Values (Taken from Tomas)
P = 0.025;       % Total phosphorus in system
b = 1.2;         % Producer maximum growth rate
q = 0.0038;      % Producer minimum P:C ratio
e_j = 0.5;       % Juvenile maximum production efficiency
e_a = 0.8;       % Adult maximum production efficiency
theta_j = 0.025; % Juvenile constant P:C ratio
theta_a = 0.03;  % Adult constant P:C ratio
delta_j = 0.06;  % Juvenile loss rate
delta_a = 0.08;  % Adult loss rate
a_h = 0.25;      % Juvenile and adult half-saturation constant
c_j = 0.5;       % Juvenile maximum ingestion rate
c_a = 0.81;      % Adult maximum ingestion rate

%% Fixed Carrying Capacity
K = 1;

%% Developmental Threshold Values
% Edit these three values to compare different developmental thresholds.
w_1 = 0.8;
w_2 = 1.4;
w_3 = 1.8;

%% Constant Non-dimensional Parameters
kappa_j = (theta_j*b)/(c_j*q);
kappa_a = (theta_a*b)/(c_a*q);
r_j = q/theta_j;
alpha = (c_j*e_a)/b;
beta = (c_a*e_a)/b;
chi_j = delta_j/b;
chi_a = delta_a/b;

%% w-dependent Non-dimensional Parameters
rho_1 = (e_j*c_j)/(w_1*b);
rho_2 = (e_j*c_j)/(w_2*b);
rho_3 = (e_j*c_j)/(w_3*b);

params_1 = {kappa_j, kappa_a, r_j, alpha, beta, ...
    chi_j, chi_a, rho_1, a_h, P, q, c_j, c_a, b};

params_2 = {kappa_j, kappa_a, r_j, alpha, beta, ...
    chi_j, chi_a, rho_2, a_h, P, q, c_j, c_a, b};

params_3 = {kappa_j, kappa_a, r_j, alpha, beta, ...
    chi_j, chi_a, rho_3, a_h, P, q, c_j, c_a, b};

%% Fixed Delay in Developmental Time
lag = 1;

%% Dimensional History
x_dim = 0.5;
J_dim = 0.125;
A_dim = 0.125;

% Each w value requires its own chronological-clock history because rho
% changes with w.
historyFcn_1 = @(phi) history( ...
    phi, params_1, K, x_dim, J_dim, A_dim);

historyFcn_2 = @(phi) history( ...
    phi, params_2, K, x_dim, J_dim, A_dim);

historyFcn_3 = @(phi) history( ...
    phi, params_3, K, x_dim, J_dim, A_dim);

%% Developmental-time Span and Solver Options
phi_span = [0, 500];
options = ddeset('RelTol', 1e-6, 'AbsTol', 1e-9);

%% Solve for the Three w Values
sol_w_1 = dde23( ...
    @(phi,y,Z) dde(phi, y, Z, params_1, K), ...
    lag, historyFcn_1, phi_span, options);

sol_w_2 = dde23( ...
    @(phi,y,Z) dde(phi, y, Z, params_2, K), ...
    lag, historyFcn_2, phi_span, options);

sol_w_3 = dde23( ...
    @(phi,y,Z) dde(phi, y, Z, params_3, K), ...
    lag, historyFcn_3, phi_span, options);

%% Non-dimensional Outputs for w_1
x_1_non_dim = sol_w_1.y(1,:);
J_1_non_dim = sol_w_1.y(2,:);
A_1_non_dim = sol_w_1.y(3,:);
t_1_non_dim = sol_w_1.y(4,:);

%% Dimensional Outputs for w_1
x_1 = x_1_non_dim*K;
J_1 = (J_1_non_dim*K*b)/c_j;
A_1 = (A_1_non_dim*K*b)/c_a;
t_1 = t_1_non_dim/b;

%% Population Plot for w_1
adult_fraction_1 = A_1./max(J_1 + A_1, eps);

figure(1);

yyaxis left
plot(t_1, J_1 + A_1, 'k-', 'LineWidth', 1.5);
hold on
plot(t_1, x_1, 'k--', 'LineWidth', 1.5);
ylabel('$\mathrm{Density\;(mg\,C\,L^{-1})}$');
ylim([0 2]);

yyaxis right
plot(t_1, adult_fraction_1, 'r-', 'LineWidth', 1.5);
ylabel('Adult proportion');
ylim([0 1]);

xlabel('$t\;(\mathrm{days})$');
title(sprintf('$K=%.3g,\\quad w=%.3g$', K, w_1));
xlim([0 300]);

legend({'Grazer $(J+A)$', 'Producer', 'Adult proportion'}, ...
    'Location', 'best');

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
box on

%% Non-dimensional Outputs for w_2
x_2_non_dim = sol_w_2.y(1,:);
J_2_non_dim = sol_w_2.y(2,:);
A_2_non_dim = sol_w_2.y(3,:);
t_2_non_dim = sol_w_2.y(4,:);

%% Dimensional Outputs for w_2
x_2 = x_2_non_dim*K;
J_2 = (J_2_non_dim*K*b)/c_j;
A_2 = (A_2_non_dim*K*b)/c_a;
t_2 = t_2_non_dim/b;

%% Population Plot for w_2
adult_fraction_2 = A_2./max(J_2 + A_2, eps);

figure(2);

yyaxis left
plot(t_2, J_2 + A_2, 'k-', 'LineWidth', 1.5);
hold on
plot(t_2, x_2, 'k--', 'LineWidth', 1.5);
ylabel('$\mathrm{Density\;(mg\,C\,L^{-1})}$');
ylim([0 2]);

yyaxis right
plot(t_2, adult_fraction_2, 'r-', 'LineWidth', 1.5);
ylabel('Adult proportion');
ylim([0 1]);

xlabel('$t\;(\mathrm{days})$');
title(sprintf('$K=%.3g,\\quad w=%.3g$', K, w_2));
xlim([0 300]);

legend({'Grazer $(J+A)$', 'Producer', 'Adult proportion'}, ...
    'Location', 'best');

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
box on

%% Non-dimensional Outputs for w_3
x_3_non_dim = sol_w_3.y(1,:);
J_3_non_dim = sol_w_3.y(2,:);
A_3_non_dim = sol_w_3.y(3,:);
t_3_non_dim = sol_w_3.y(4,:);

%% Dimensional Outputs for w_3
x_3 = x_3_non_dim*K;
J_3 = (J_3_non_dim*K*b)/c_j;
A_3 = (A_3_non_dim*K*b)/c_a;
t_3 = t_3_non_dim/b;

%% Population Plot for w_3
adult_fraction_3 = A_3./max(J_3 + A_3, eps);

figure(3);

yyaxis left
plot(t_3, J_3 + A_3, 'k-', 'LineWidth', 1.5);
hold on
plot(t_3, x_3, 'k--', 'LineWidth', 1.5);
ylabel('$\mathrm{Density\;(mg\,C\,L^{-1})}$');
ylim([0 2]);

yyaxis right
plot(t_3, adult_fraction_3, 'r-', 'LineWidth', 1.5);
ylabel('Adult proportion');
ylim([0 1]);

xlabel('$t\;(\mathrm{days})$');
title(sprintf('$K=%.3g,\\quad w=%.3g$', K, w_3));
xlim([0 300]);

legend({'Grazer $(J+A)$', 'Producer', 'Adult proportion'}, ...
    'Location', 'best');

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
box on

%% Plot Q over Time for the Three w Values
Q_1 = (P - theta_j.*J_1 - theta_a.*A_1)./x_1;
Q_2 = (P - theta_j.*J_2 - theta_a.*A_2)./x_2;
Q_3 = (P - theta_j.*J_3 - theta_a.*A_3)./x_3;

figure(4);
clf

semilogy(t_1, Q_1, 'k:', 'LineWidth', 1.5);
hold on
semilogy(t_2, Q_2, 'k--', 'LineWidth', 1.5);
semilogy(t_3, Q_3, 'k-', 'LineWidth', 1.5);

xlabel('$t\;(\mathrm{days})$');
ylabel('$Q\;(\mathrm{mg\,P\,mg^{-1}\,C})$');
title(sprintf('$K=%.3g$', K));

xlim([0 300]);
ylim([5e-3 2e1]);

legend({ ...
    sprintf('$w=%.3g$', w_1), ...
    sprintf('$w=%.3g$', w_2), ...
    sprintf('$w=%.3g$', w_3)}, ...
    'Location', 'northeast');

box on

%% Stage-timing Information for w_1
[t_1_days, tau_1_days, estimated_period_1] = ...
    stage_timing_information( ...
        sol_w_1, t_1_non_dim, J_1, A_1, b, 100, 'w_1');

%% Stage-timing Plot for w_1
figure(5);
clf

plot(t_1_days, tau_1_days, 'k-', 'LineWidth', 1.5);
hold on

if ~isnan(estimated_period_1)
    yline(estimated_period_1, '--', ...
        'Color', [0.65 0.65 0.65], ...
        'LineWidth', 1.5);

    legend({'$\tau(t)$', '$\mathrm{Estimated\ cycle\ period}$'}, ...
        'Location', 'best');
else
    legend({'$\tau(t)$'}, 'Location', 'best');
end

xlabel('$t\;(\mathrm{days})$');
ylabel('Stage duration and cycle period $(\mathrm{days})$');
title(sprintf('$K=%.3g,\\quad w=%.3g$', K, w_1));
xlim([0 1000]);

upper_limit_1 = 1.1*max( ...
    [tau_1_days, estimated_period_1], [], 'omitnan');
ylim([0 upper_limit_1]);
box on

fprintf('\nw_1 = %.3f at K = %.3f\n', w_1, K);
fprintf('Estimated cycle period: %.3f days\n', estimated_period_1);
fprintf('Stage-duration range: %.3f to %.3f days\n', ...
    min(tau_1_days), max(tau_1_days));

%% Stage-timing Information for w_2
[t_2_days, tau_2_days, estimated_period_2] = ...
    stage_timing_information( ...
        sol_w_2, t_2_non_dim, J_2, A_2, b, 100, 'w_2');

%% Stage-timing Plot for w_2
figure(6);
clf

plot(t_2_days, tau_2_days, 'k-', 'LineWidth', 1.5);
hold on

if ~isnan(estimated_period_2)
    yline(estimated_period_2, '--', ...
        'Color', [0.65 0.65 0.65], ...
        'LineWidth', 1.5);

    legend({'$\tau(t)$', '$\mathrm{Estimated\ cycle\ period}$'}, ...
        'Location', 'best');
else
    legend({'$\tau(t)$'}, 'Location', 'best');
end

xlabel('$t\;(\mathrm{days})$');
ylabel('Stage duration and cycle period $(\mathrm{days})$');
title(sprintf('$K=%.3g,\\quad w=%.3g$', K, w_2));
xlim([0 1000]);

upper_limit_2 = 1.1*max( ...
    [tau_2_days, estimated_period_2], [], 'omitnan');
ylim([0 upper_limit_2]);
box on

fprintf('\nw_2 = %.3f at K = %.3f\n', w_2, K);
fprintf('Estimated cycle period: %.3f days\n', estimated_period_2);
fprintf('Stage-duration range: %.3f to %.3f days\n', ...
    min(tau_2_days), max(tau_2_days));

%% Stage-timing Information for w_3
[t_3_days, tau_3_days, estimated_period_3] = ...
    stage_timing_information( ...
        sol_w_3, t_3_non_dim, J_3, A_3, b, 100, 'w_3');

%% Stage-timing Plot for w_3
figure(7);
clf

plot(t_3_days, tau_3_days, 'k-', 'LineWidth', 1.5);
hold on

if ~isnan(estimated_period_3)
    yline(estimated_period_3, '--', ...
        'Color', [0.65 0.65 0.65], ...
        'LineWidth', 1.5);

    legend({'$\tau(t)$', '$\mathrm{Estimated\ cycle\ period}$'}, ...
        'Location', 'best');
else
    legend({'$\tau(t)$'}, 'Location', 'best');
end

xlabel('$t\;(\mathrm{days})$');
ylabel('Stage duration and cycle period $(\mathrm{days})$');
title(sprintf('$K=%.3g,\\quad w=%.3g$', K, w_3));
xlim([0 1000]);

upper_limit_3 = 1.1*max( ...
    [tau_3_days, estimated_period_3], [], 'omitnan');
ylim([0 upper_limit_3]);
box on

fprintf('\nw_3 = %.3f at K = %.3f\n', w_3, K);
fprintf('Estimated cycle period: %.3f days\n', estimated_period_3);
fprintf('Stage-duration range: %.3f to %.3f days\n', ...
    min(tau_3_days), max(tau_3_days));

%% Model Functions

function output = dde(~, y, Z, params, K)
    % Unpack parameters
    [kappa_j, kappa_a, r_j, alpha, beta, ...
        chi_j, chi_a, rho, a_h, P, q, ~, ~, ~] = ...
        deal(params{:});

    lambda = a_h/K;
    Pi = P/(q*K);

    % Preallocate derivative vector
    output = zeros(size(y));

    % Extract variables
    x = y(1);
    J = y(2);
    A = y(3);
    T = y(4);

    % Define Q, h_j, and D
    Q = (Pi - kappa_j*J - kappa_a*A)/x;
    h_j = tanh(r_j*Q);
    D = h_j*(x/(lambda + x));

    % Survival probability
    S_phi = exp(-chi_j*(T - Z(4,1)));

    % Population equations
    output(1) = ...
        (x/(rho*D))*(1-x)*(1-(1/Q)) ...
        - (J+A)/(rho*h_j);

    output(2) = ...
        (alpha/rho)*(A - Z(3,1)*S_phi) ...
        - (chi_j/(rho*D))*J;

    output(3) = ...
        (beta/rho)*Z(3,1)*S_phi ...
        - (chi_a/(rho*D))*A;

    % Chronological-clock equation
    output(4) = 1/(rho*D);
end


function yhist = history(phi, params, K, x_dim, J_dim, A_dim)
    % Unpack parameters needed for the history
    [kappa_j, kappa_a, r_j, ~, ~, ~, ~, rho, ...
        a_h, P, q, c_j, c_a, b] = deal(params{:});

    lambda = a_h/K;
    Pi = P/(q*K);

    % Construct constant biological histories
    x0 = x_dim/K;
    J0 = (J_dim*c_j)/(b*K);
    A0 = (A_dim*c_a)/(b*K);

    % Evaluate the constant-history development rate
    Q0 = (Pi - kappa_j*J0 - kappa_a*A0)/x0;
    h_j0 = tanh(r_j*Q0);
    F0 = x0/(lambda + x0);
    D0 = h_j0*F0;

    % Check that the developmental clock is well defined
    if x0 <= 0 || Q0 <= 0 || D0 <= 0
        error('The selected history is not biologically feasible.');
    end

    nPhi = numel(phi);

    % Constant biological histories and linear clock history
    yhist = [
        x0*ones(1,nPhi);
        J0*ones(1,nPhi);
        A0*ones(1,nPhi);
        phi/(rho*D0)
    ];
end


function [t_days, tau_days, estimated_period] = ...
    stage_timing_information( ...
        sol, T_non_dim, J_dim, A_dim, b, transient_end, case_name)

    phi = sol.x;

    % T(phi-1) is available from the numerical solution for phi >= 1.
    valid = phi >= 1;
    phi_valid = phi(valid);

    T_now = T_non_dim(valid);
    T_delayed = deval(sol, phi_valid - 1, 4);

    % Convert nondimensional chronological time to days.
    t_days = T_now/b;
    tau_days = (T_now - T_delayed)/b;

    % Estimate the late-time cycle period from total grazer biomass.
    grazer = J_dim(valid) + A_dim(valid);
    late = t_days >= transient_end;

    t_late = t_days(late);
    grazer_late = grazer(late);

    if numel(t_late) < 3 || range(grazer_late) <= eps
        estimated_period = NaN;
        warning('Not enough late-time variation for %s.', case_name);
        return
    end

    minimum_prominence = 0.05*range(grazer_late);

    [~, peak_times] = findpeaks( ...
        grazer_late, t_late, ...
        'MinPeakProminence', minimum_prominence);

    if numel(peak_times) >= 2
        estimated_period = median(diff(peak_times));
    else
        estimated_period = NaN;
        warning( ...
            'Not enough late-time peaks to estimate the period for %s.', ...
            case_name);
    end
end
