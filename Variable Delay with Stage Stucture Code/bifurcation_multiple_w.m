%% Basic Bifurcation Analysis of Combined Model
%% Varying values of W for a bifurcation analysis
clear;
close all;
clc;

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

%% Parameter values (taken from Tomas)
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

%% Stage-duration parameters
% Specify the three developmental thresholds here.
w_values = [1.00, 1.23, 1.46];
n_w = numel(w_values);

%% Constant nondimensional parameters that do not depend on w
kappa_j = (theta_j*b)/(c_j*q);
kappa_a = (theta_a*b)/(c_a*q);
r_j = q/theta_j;
alpha = (c_j*e_a)/b;
beta = (c_a*e_a)/b;
chi_j = delta_j/b;
chi_a = delta_a/b;

%% Delay
lag = 1;

%% History
x_dim = 0.3125;
J_dim = 0.09;
A_dim = 0.01;

%% K range
K_range = 0.1:0.05:3;

%% Developmental-time span
phi_span = [0, 5000];

%% Solver options
options = ddeset('RelTol', 1e-6, 'AbsTol', 1e-9);

%% Late-time and period-estimation settings
late_window_days = 5000;
sample_step_days = 0.25;

% A period is reported only if the final relative amplitude exceeds this
% threshold. This prevents numerical noise near an equilibrium from being
% interpreted as a cycle.
cycle_relative_amplitude_threshold = 1e-3;

%% Storage
x_min_storage = nan(n_w, numel(K_range));
x_max_storage = nan(n_w, numel(K_range));

grazer_min_storage = nan(n_w, numel(K_range));
grazer_max_storage = nan(n_w, numel(K_range));

period_storage = nan(n_w, numel(K_range));
relative_amplitude_storage = nan(n_w, numel(K_range));
number_of_peaks_storage = zeros(n_w, numel(K_range));

%% Bifurcation loops
for w_index = 1:n_w

    w = w_values(w_index);
    rho = (e_j*c_j)/(w*b);

    params = { ...
        kappa_j, kappa_a, r_j, alpha, beta, ...
        chi_j, chi_a, rho, a_h, P, q, c_j, c_a, b};

    fprintf('\nStarting scan for w = %.3f\n', w);

    for K_index = 1:numel(K_range)

        K = K_range(K_index);

        fprintf('  Current K value: %.3f\n', K);

        historyFcn = @(phi) history( ...
            phi, params, K, x_dim, J_dim, A_dim);

        sol = dde23( ...
            @(phi,y,Z) dde(phi, y, Z, params, K), ...
            lag, historyFcn, phi_span, options);

        %% Convert to dimensional variables
        t_days = sol.y(4,:)/b;

        x_dimensional = K*sol.y(1,:);
        J_dimensional = (b*K/c_j)*sol.y(2,:);
        A_dimensional = (b*K/c_a)*sol.y(3,:);

        grazer_dimensional = J_dimensional + A_dimensional;

        %% Remove duplicate chronological times, if present
        [t_days, unique_indices] = unique(t_days, 'stable');
        x_dimensional = x_dimensional(unique_indices);
        grazer_dimensional = grazer_dimensional(unique_indices);

        finite_values = ...
            isfinite(t_days) ...
            & isfinite(x_dimensional) ...
            & isfinite(grazer_dimensional);

        t_days = t_days(finite_values);
        x_dimensional = x_dimensional(finite_values);
        grazer_dimensional = grazer_dimensional(finite_values);

        if numel(t_days) < 10 || any(diff(t_days) <= 0)
            warning(['Insufficient chronological-time output at ', ...
                'w = %.3f, K = %.3f.'], w, K);
            continue
        end

        %% Retain only the final chronological-time window
        late_start = max(t_days(end)-late_window_days, t_days(1));
        late = t_days >= late_start;

        if nnz(late) < 10
            warning('Too few late-time points at w = %.3f, K = %.3f.', ...
                w, K);
            continue
        end

        %% Store late-time extrema
        x_min_storage(w_index,K_index) = ...
            min(x_dimensional(late));
        x_max_storage(w_index,K_index) = ...
            max(x_dimensional(late));

        grazer_min_storage(w_index,K_index) = ...
            min(grazer_dimensional(late));
        grazer_max_storage(w_index,K_index) = ...
            max(grazer_dimensional(late));

        %% Estimate the cycle period in chronological days
        %
        % dde23 returns a nonuniform mesh. Interpolate onto a uniform grid
        % before identifying peaks.
        t_uniform = late_start:sample_step_days:t_days(end);

        if numel(t_uniform) < 10
            warning(['Uniform late-time grid is too short at ', ...
                'w = %.3f, K = %.3f.'], w, K);
            continue
        end

        x_uniform = interp1( ...
            t_days, x_dimensional, t_uniform, 'pchip');

        grazer_uniform = interp1( ...
            t_days, grazer_dimensional, t_uniform, 'pchip');

        resource_amplitude = max(x_uniform)-min(x_uniform);
        grazer_amplitude = max(grazer_uniform)-min(grazer_uniform);

        relative_resource_amplitude = ...
            resource_amplitude/max(abs(mean(x_uniform)), 1e-12);

        relative_grazer_amplitude = ...
            grazer_amplitude/max(abs(mean(grazer_uniform)), 1e-12);

        relative_amplitude = max( ...
            relative_resource_amplitude, ...
            relative_grazer_amplitude);

        relative_amplitude_storage(w_index,K_index) = ...
            relative_amplitude;

        if relative_amplitude >= cycle_relative_amplitude_threshold
            [estimated_period, number_of_peaks] = ...
                estimate_cycle_period(t_uniform, x_uniform);

            period_storage(w_index,K_index) = estimated_period;
            number_of_peaks_storage(w_index,K_index) = ...
                number_of_peaks;
        end

        if isfinite(period_storage(w_index,K_index))
            fprintf('    Estimated period: %.3f days using %d peaks\n', ...
                period_storage(w_index,K_index), ...
                number_of_peaks_storage(w_index,K_index));
        else
            fprintf('    No resolved late-time cycle\n');
        end
    end
end

%% Plot bifurcation diagrams and cycle period
colors = lines(n_w);
legend_labels = arrayfun( ...
    @(value) sprintf('$w=%.2f$', value), ...
    w_values, 'UniformOutput', false);

figure('Color', 'w');
tiledlayout(1, 3, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

nexttile
hold on

resource_legend_handles = gobjects(n_w,1);

for w_index = 1:n_w
    resource_legend_handles(w_index) = plot( ...
        K_range, x_max_storage(w_index,:), '.', ...
        'Color', colors(w_index,:), ...
        'MarkerSize', 10);

    plot(K_range, x_min_storage(w_index,:), '.', ...
        'Color', colors(w_index,:), ...
        'MarkerSize', 10, ...
        'HandleVisibility', 'off');
end

xlabel('$K\;(\mathrm{mg\,C\,L^{-1}})$');
ylabel('$\mathrm{Resource\ density\;(mg\,C\,L^{-1})}$');
title('Resource bifurcation diagram');
legend(resource_legend_handles, legend_labels, ...
    'Location', 'best');
grid on
box on

nexttile
hold on

grazer_legend_handles = gobjects(n_w,1);

for w_index = 1:n_w
    grazer_legend_handles(w_index) = plot( ...
        K_range, grazer_max_storage(w_index,:), '.', ...
        'Color', colors(w_index,:), ...
        'MarkerSize', 10);

    plot(K_range, grazer_min_storage(w_index,:), '.', ...
        'Color', colors(w_index,:), ...
        'MarkerSize', 10, ...
        'HandleVisibility', 'off');
end

xlabel('$K\;(\mathrm{mg\,C\,L^{-1}})$');
ylabel('$J+A\;(\mathrm{mg\,C\,L^{-1}})$');
title('Total grazer bifurcation diagram');
legend(grazer_legend_handles, legend_labels, ...
    'Location', 'best');
grid on
box on

nexttile
hold on

period_legend_handles = gobjects(n_w,1);

for w_index = 1:n_w
    period_legend_handles(w_index) = plot( ...
        K_range, period_storage(w_index,:), '.-', ...
        'Color', colors(w_index,:), ...
        'MarkerSize', 12, ...
        'LineWidth', 1.2);
end

xlabel('$K\;(\mathrm{mg\,C\,L^{-1}})$');
ylabel('$\mathrm{Estimated\ cycle\ period\;(days)}$');
title('Cycle period');
legend(period_legend_handles, legend_labels, ...
    'Location', 'best');
grid on
box on

%% Print maximum resolved period for each w
for w_index = 1:n_w
    finite_periods = isfinite(period_storage(w_index,:));

    if any(finite_periods)
        finite_indices = find(finite_periods);
        [maximum_period, local_index] = ...
            max(period_storage(w_index,finite_periods));
        maximum_index = finite_indices(local_index);

        fprintf('\nw = %.3f\n', w_values(w_index));
        fprintf('  Maximum estimated period: %.6f days\n', ...
            maximum_period);
        fprintf('  K at maximum period: %.6f\n', ...
            K_range(maximum_index));
    else
        fprintf('\nw = %.3f: no cycles were resolved.\n', ...
            w_values(w_index));
    end
end

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
    Q = (Pi-kappa_j*J-kappa_a*A)/x;
    h_j = tanh(r_j*Q);
    D = h_j*(x/(lambda+x));

    % Survival probability
    S_phi = exp(-chi_j*(T-Z(4,1)));

    % Population equations
    output(1) = ...
        (x/(rho*D))*(1-x)*(1-(1/Q)) ...
        -(J+A)/(rho*h_j);

    output(2) = ...
        (alpha/rho)*(A-Z(3,1)*S_phi) ...
        -(chi_j/(rho*D))*J;

    output(3) = ...
        (beta/rho)*Z(3,1)*S_phi ...
        -(chi_a/(rho*D))*A;

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
    Q0 = (Pi-kappa_j*J0-kappa_a*A0)/x0;
    h_j0 = tanh(r_j*Q0);
    F0 = x0/(lambda+x0);
    D0 = h_j0*F0;

    % Check that the developmental clock is well defined
    if x0 <= 0 || Q0 <= 0 || D0 <= 0
        error('The selected history is not biologically feasible.');
    end

    nPhi = numel(phi);

    % Constant biological histories and linear clock history
    yhist = [ ...
        x0*ones(1,nPhi); ...
        J0*ones(1,nPhi); ...
        A0*ones(1,nPhi); ...
        phi/(rho*D0)];
end

function [period, number_of_peaks] = ...
        estimate_cycle_period(time, signal)

    period = NaN;
    number_of_peaks = 0;

    if numel(signal) < 5
        return
    end

    % Identify local maxima without requiring the Signal Processing Toolbox.
    peak_indices = find( ...
        signal(2:end-1) > signal(1:end-2) ...
        & signal(2:end-1) >= signal(3:end))+1;

    signal_amplitude = max(signal)-min(signal);

    if signal_amplitude <= 0
        return
    end

    % Retain peaks that lie in the upper quarter of the oscillation range.
    peak_threshold = min(signal)+0.75*signal_amplitude;
    peak_indices = peak_indices( ...
        signal(peak_indices) >= peak_threshold);

    number_of_peaks = numel(peak_indices);

    if number_of_peaks >= 3
        period = median(diff(time(peak_indices)));
    end
end
