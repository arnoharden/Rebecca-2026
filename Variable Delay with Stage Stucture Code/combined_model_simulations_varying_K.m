%% Creating Tomas-esque Figures Using Combined Stoichiometric Model

%Daphnia = Black
%Algea = Gray
%Cycle period = dashed line
%Stage duration = solid line

clear;
close all;
clc;

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

%%  Parameter Values (Taken from Tomas)
P = 0.025;       %Total phosphorous in system
b = 1.2;         %Producer max growth rate
q = 0.0038;      %Producer min P:X
e_j = 0.5;       %Juevenile max production efficiency
e_a = 0.8;       %Adult max production efficiency
theta_j = 0.025; %Juvenile constant P:C ratio
theta_a = 0.03;  %Adult constant P:C ratio
delta_j = 0.06;  %Juvenile loss rate
delta_a = 0.08;  %Adult loss rate
a_h = 0.25;      %Juvenile and adult half-saturation constant
c_j = 0.5;       %Juvenile max ingestion rate
c_a = 0.81;      %Adult max ingestion rate 

%% Stage Duration Parameter
w = 1;

%% Constant Non-dimensional parameters
kappa_j = (theta_j*b)/(c_j*q);
kappa_a = (theta_a*b)/(c_a*q);
r_j = q/theta_j;
alpha = (c_j*e_a)/b;
beta = (c_a*e_a)/b;
chi_j = delta_j/b;
chi_a = delta_a/b;
rho = (e_j*c_j)/(w*b);

params = {kappa_j, kappa_a, r_j, alpha, beta, chi_j, chi_a, rho, a_h, P, q, c_j, c_a, b};

%% Carrying Capacity Variables
K_1 = 0.5;
K_2 = 1;
K_3 = 2;

%% Delays
lag = 1;

%% History
x_dim = 0.5;
J_dim = 0.125;
A_dim = 0.125;

historyFcn_1 = @(phi) history(phi, params, K_1, x_dim, J_dim, A_dim);
historyFcn_2 = @(phi) history(phi, params, K_2, x_dim, J_dim, A_dim);
historyFcn_3 = @(phi) history(phi, params, K_3, x_dim, J_dim, A_dim);

%% Phi tspans
phi_span = [0,500];

%% options
options = ddeset('RelTol', 1e-6, 'AbsTol', 1e-9);

%% Solving for varying K values
sol_K_1 = dde23(@(phi,y,Z) dde(phi, y, Z, params, K_1), lag, historyFcn_1, phi_span, options);
sol_K_2 = dde23(@(phi,y,Z) dde(phi, y, Z, params, K_2), lag, historyFcn_2, phi_span, options);
sol_K_3 = dde23(@(phi,y,Z) dde(phi, y, Z, params, K_3), lag, historyFcn_3, phi_span, options);

%% Non-dimensional outputs for K = 0.5
x_1_non_dim = sol_K_1.y(1,:);
J_1_non_dim = sol_K_1.y(2,:);
A_1_non_dim = sol_K_1.y(3,:);
t_1_non_dim = sol_K_1.y(4,:);

%% Dimensional outputs for K = 0.5
x_1 = x_1_non_dim*K_1;
J_1 = (J_1_non_dim*K_1*b)/c_j;
A_1 = (A_1_non_dim*K_1*b)/c_a;
t_1 = t_1_non_dim/b;

%% Plotting for K = 0.5
adult_fraction_1 = A_1./max(J_1 + A_1, eps);

figure(1);

yyaxis left
plot(t_1, J_1 + A_1, 'k-', 'LineWidth', 1.5);
hold on
plot(t_1, x_1, '--k', 'LineWidth', 1.5);
ylabel('$\mathrm{Density\;(mg\,C\,L^{-1})}$');
ylim([0 2]);

yyaxis right
plot(t_1, adult_fraction_1, 'r-', 'LineWidth', 1.5);
ylabel('Adult Proportion');
ylim([0 1]);

xlabel('$t\;(\mathrm{days})$');
title('$(a)\quad K=0.5$');
xlim([0 300]);

legend({'Grazer (J+A)', 'Producer', 'Adult Proportion'}, ...
    'Location', 'best');

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';

%% Non-dimensional outputs for K = 1
x_2_non_dim = sol_K_2.y(1,:);
J_2_non_dim = sol_K_2.y(2,:);
A_2_non_dim = sol_K_2.y(3,:);
t_2_non_dim = sol_K_2.y(4,:);

%% Dimensional outputs for K = 1
x_2 = x_2_non_dim*K_2;
J_2 = (J_2_non_dim*K_2*b)/c_j;
A_2 = (A_2_non_dim*K_2*b)/c_a;
t_2 = t_2_non_dim/b;

%% Plotting for K = 1
adult_fraction_2 = A_2./max(J_2 + A_2, eps);

figure(2);

yyaxis left
plot(t_2, J_2 + A_2, 'k-', 'LineWidth', 1.5);
hold on
plot(t_2, x_2, '--k', 'LineWidth', 1.5);
ylabel('$\mathrm{Density\;(mg\,C\,L^{-1})}$');
ylim([0 2]);

yyaxis right
plot(t_2, adult_fraction_2, 'r-', 'LineWidth', 1.5);
ylabel('Adult Proportion');
ylim([0 1]);

xlabel('$t\;(\mathrm{days})$');
title('$(b)\quad K=1$');
xlim([0 300]);

legend({'Grazer (J+A)', 'Producer', 'Adult Proportion'}, ...
    'Location', 'best');

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';

%% Non-dimensional outputs for K = 2
x_3_non_dim = sol_K_3.y(1,:);
J_3_non_dim = sol_K_3.y(2,:);
A_3_non_dim = sol_K_3.y(3,:);
t_3_non_dim = sol_K_3.y(4,:);

%% Dimensional outputs for K = 2
x_3 = x_3_non_dim*K_3;
J_3 = (J_3_non_dim*K_3*b)/c_j;
A_3 = (A_3_non_dim*K_3*b)/c_a;
t_3 = t_3_non_dim/b;

%% Plotting for K = 2
adult_fraction_3 = A_3./max(J_3 + A_3, eps);

figure(3);

yyaxis left
plot(t_3, J_3 + A_3, 'k-', 'LineWidth', 1.5);
hold on
plot(t_3, x_3, '--k', 'LineWidth', 1.5);
ylabel('$\mathrm{Density\;(mg\,C\,L^{-1})}$');
ylim([0 2]);

yyaxis right
plot(t_3, adult_fraction_3, 'r-', 'LineWidth', 1.5);
ylabel('Adult Proportion');
ylim([0 1]);

xlabel('$t\;(\mathrm{days})$');
title('$(c)\quad K=2$');
xlim([0 300]);

legend({'Grazer (J+A)', 'Producer', 'Adult Proportion'}, ...
    'Location', 'best');

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';

%% Plotting Q over time

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

xlim([0 300]);
ylim([5e-3 2e1]);

legend({ ...
    '$K=0.5\;\mathrm{mg\,C\,L^{-1}}$', ...
    '$K=1\;\mathrm{mg\,C\,L^{-1}}$', ...
    '$K=2\;\mathrm{mg\,C\,L^{-1}}$'}, ...
    'Location', 'northeast');

box on

%% Variable delay for K = 0.5

phi_1 = sol_K_1.x;

% T(phi-1) is available for phi >= 1
valid_1 = phi_1 >= 1;
phi_valid_1 = phi_1(valid_1);

T_1_now = t_1_non_dim(valid_1);
T_1_delayed = deval(sol_K_1, phi_valid_1 - 1, 4);

% Convert to chronological days
t_1_days = T_1_now/b;
tau_1_days = (T_1_now - T_1_delayed)/b;

%% Estimate the late-time cycle period

grazer_1 = J_1(valid_1) + A_1(valid_1);

transient_end = 100;
late_1 = t_1_days >= transient_end;

t_1_late = t_1_days(late_1);
grazer_1_late = grazer_1(late_1);

minimum_prominence_1 = 0.05*range(grazer_1_late);

[~, peak_times_1] = findpeaks( ...
    grazer_1_late, ...
    t_1_late, ...
    'MinPeakProminence', minimum_prominence_1);

if numel(peak_times_1) >= 2
    estimated_period_1 = median(diff(peak_times_1));
else
    estimated_period_1 = NaN;
    warning('Not enough late-time peaks for K_1.');
end

%% Stage-timing plot for K_1

figure(5);
clf

plot(t_1_days, tau_1_days, ...
    'k-', 'LineWidth', 1.5);
hold on

if ~isnan(estimated_period_1)
    yline(estimated_period_1, ...
        '--', ...
        'Color', [0.65 0.65 0.65], ...
        'LineWidth', 1.5);

    legend({'$\tau(t)$', '$\mathrm{Estimated\ cycle\ period}$'}, ...
        'Location', 'best');
else
    legend({'$\tau(t)$'}, 'Location', 'best');
end

xlabel('$t\;(\mathrm{days})$');
ylabel('Stage duration and cycle period $(\mathrm{days})$');
title(sprintf('$K=%.3g$', K_1));

xlim([0 500]);

upper_limit_1 = 1.1*max( ...
    [tau_1_days, estimated_period_1], [], 'omitnan');
ylim([0 upper_limit_1]);

box on

fprintf('\nK_1 = %.3f\n', K_1);
fprintf('Estimated cycle period: %.3f days\n', ...
    estimated_period_1);
fprintf('Stage-duration range: %.3f to %.3f days\n', ...
    min(tau_1_days), max(tau_1_days));


%% Plotting variable delay for K = 1

phi_2 = sol_K_2.x;

% T(phi-1) is available for phi >= 1
valid = phi_2 >= 1;
phi_valid_2 = phi_2(valid);

T_2_now = t_2_non_dim(valid);
T_2_delayed = deval(sol_K_2, phi_valid_2 - 1, 4);

% Convert to chronological days
t_2_days = T_2_now/b;
tau_2_days = (T_2_now - T_2_delayed)/b;

%% Estimate the late-time cycle period

grazer_2 = J_2(valid) + A_2(valid);

transient_end = 100;
late = t_2_days >= transient_end;

t_late = t_2_days(late);
grazer_late = grazer_2(late);

minimum_prominence = 0.05*range(grazer_late);

[~, peak_times] = findpeaks( ...
    grazer_late, ...
    t_late, ...
    'MinPeakProminence', minimum_prominence);

if numel(peak_times) >= 2
    estimated_period = median(diff(peak_times));
else
    estimated_period = NaN;
    warning('Not enough late-time peaks to estimate the cycle period.');
end

%% Stage-timing plot

figure(6);
clf

plot(t_2_days, tau_2_days, ...
    'k-', 'LineWidth', 1.5);
hold on

if ~isnan(estimated_period)
    yline(estimated_period, ...
        '--', ...
        'Color', [0.65 0.65 0.65], ...
        'LineWidth', 1.5);
end

xlabel('$t\;(\mathrm{days})$');
ylabel('Stage duration and cycle period $(\mathrm{days})$');
title(sprintf('$K=%.3g$', K_2));

legend({'$\tau(t)$', '$\mathrm{Estimated\ cycle\ period}$'}, ...
    'Location', 'best');

xlim([0 500]);

upper_limit = 1.1*max([tau_2_days, estimated_period], ...
    [], 'omitnan');
ylim([0 upper_limit]);

box on

fprintf('Estimated cycle period: %.3f days\n', estimated_period);
fprintf('Stage-duration range: %.3f to %.3f days\n', ...
    min(tau_2_days), max(tau_2_days));

%% Variable delay for K = K_3

phi_3 = sol_K_3.x;

% T(phi-1) is available for phi >= 1
valid_3 = phi_3 >= 1;
phi_valid_3 = phi_3(valid_3);

T_3_now = t_3_non_dim(valid_3);
T_3_delayed = deval(sol_K_3, phi_valid_3 - 1, 4);

% Convert to chronological days
t_3_days = T_3_now/b;
tau_3_days = (T_3_now - T_3_delayed)/b;

%% Estimate the late-time cycle period

grazer_3 = J_3(valid_3) + A_3(valid_3);

transient_end = 100;
late_3 = t_3_days >= transient_end;

t_3_late = t_3_days(late_3);
grazer_3_late = grazer_3(late_3);

minimum_prominence_3 = 0.05*range(grazer_3_late);

[~, peak_times_3] = findpeaks( ...
    grazer_3_late, ...
    t_3_late, ...
    'MinPeakProminence', minimum_prominence_3);

if numel(peak_times_3) >= 2
    estimated_period_3 = median(diff(peak_times_3));
else
    estimated_period_3 = NaN;
    warning('Not enough late-time peaks for K_3.');
end

%% Stage-timing plot for K_3

figure(7);
clf

plot(t_3_days, tau_3_days, ...
    'k-', 'LineWidth', 1.5);
hold on

if ~isnan(estimated_period_3)
    yline(estimated_period_3, ...
        '--', ...
        'Color', [0.65 0.65 0.65], ...
        'LineWidth', 1.5);

    legend({'$\tau(t)$', '$\mathrm{Estimated\ cycle\ period}$'}, ...
        'Location', 'best');
else
    legend({'$\tau(t)$'}, 'Location', 'best');
end

xlabel('$t\;(\mathrm{days})$');
ylabel('Stage duration and cycle period $(\mathrm{days})$');
title(sprintf('$K=%.3g$', K_3));

xlim([0 500]);

upper_limit_3 = 1.1*max( ...
    [tau_3_days, estimated_period_3], [], 'omitnan');
ylim([0 upper_limit_3]);

box on

fprintf('\nK_3 = %.3f\n', K_3);
fprintf('Estimated cycle period: %.3f days\n', ...
    estimated_period_3);
fprintf('Stage-duration range: %.3f to %.3f days\n', ...
    min(tau_3_days), max(tau_3_days));


function output = dde(phi, y, Z, params, K)
    %unpack params
    [kappa_j, kappa_a, r_j, alpha, beta, chi_j, chi_a, rho, a_h, P, q, ~, ~, ~] = deal(params{:});
    lambda = a_h/K;
    Pi = P/(q*K);

    %preallocate derivative vector 
    output = zeros(size(y));

    %extract variables
    x = y(1);
    J = y(2);
    A = y(3);
    t = y(4);

    %define Q, h_j, D
    Q = (Pi - kappa_j*J - kappa_a*A)/x;
    h_j = tanh(r_j * Q);
    D = h_j*(x/(lambda + x));

    %survival probability
    S_phi = exp(-chi_j*(t - Z(4,1)));

    %population equations
    output(1) = (x/(rho*D))*(1-x)*(1-(1/Q)) - (J+A)/(rho*h_j);
    output(2) = (alpha/rho)*(A - Z(3,1)*S_phi)- (chi_j/(rho*D))*J;
    output(3) = (beta/rho)*Z(3,1)*S_phi - (chi_a/(rho*D))*A;

    %time equation
    output(4) = 1/(rho*D);
end


function yhist = history(phi, params, K, x_dim, J_dim, A_dim)

    % Unpack parameters needed for the history
    [kappa_j, kappa_a, r_j, ~, ~, ~, ~, rho, ...
        a_h, P, q, c_j, c_a, b] = deal(params{:});

    lambda = a_h / K;
    Pi = P / (q * K);
    
    % Construct constant histories
    x0 = x_dim/K;
    J0 = (J_dim*c_j)/(b*K);
    A0 = (A_dim*c_a)/(b*K);

    % Evaluate the constant-history development rate
    Q0 = (Pi - kappa_j*J0 - kappa_a*A0) / x0;
    h_j0 = tanh(r_j * Q0);
    F0 = x0 / (lambda + x0);
    D0 = h_j0 * F0;

    % Check that the developmental clock is well defined
    if x0 <= 0 || Q0 <= 0 || D0 <= 0
        error('The selected history is not biologically feasible.');
    end

    nPhi = numel(phi);

    yhist = [
        x0 * ones(1,nPhi);
        J0 * ones(1,nPhi);
        A0 * ones(1,nPhi);
        phi / (rho*D0)
    ];
end
