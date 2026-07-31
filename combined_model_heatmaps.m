%% Heatmaps for the Combined Stoichiometric Variable-Delay Model
% This script saves a mat file witht the neccesary data to constuct the
% heatmaps. Use plot_combined_model_heatmaps.m for plotting
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

%% Constant nondimensional parameters
kappa_j = (theta_j*b)/(c_j*q);
kappa_a = (theta_a*b)/(c_a*q);
r_j = q/theta_j;
alpha = (c_j*e_a)/b;
beta = (c_a*e_a)/b;
chi_j = delta_j/b;
chi_a = delta_a/b;

params = { ...
    kappa_j, kappa_a, r_j, alpha, beta, ...
    chi_j, chi_a, a_h, P, q, c_j, c_a, e_j, b};

%% K and w parameter ranges
% K=0 and w=0 are excluded because both occur in denominators.
% Use a coarser grid for testing and a finer grid for the final figure.
K_range = 0.1:0.1:2;
w_range = 0.1:0.1:2;

nK = numel(K_range);
nW = numel(w_range);

%% Empty heatmap matrices
% Rows correspond to w and columns correspond to K.
grazer_amplitude = nan(nW,nK);
cycle_period_days = nan(nW,nK);
mean_stage_duration_days = nan(nW,nK);

% Classification codes:
% 0 = solver failure or infeasible history
% 1 = extinction
% 2 = non-cyclic persistence
% 3 = limit cycle
outcome_code = zeros(nW,nK);

%% Delay
lag = 1;

%% Dimensional constant history
x_dim = 0.3125;
J_dim = 0.125;
A_dim = 0.125;

%% Developmental-time span and solver options
phi_span = [0,5000];

options = ddeset( ...
    'RelTol',1e-7, ...
    'AbsTol',1e-10);

%% Analysis settings
late_window_days = 5000;
sample_step_days = 0.25;

% A solution is treated as extinct when its late-time total grazer
% biomass remains below this dimensional threshold.
extinction_threshold = 1e-6;

% A cycle must have a non-negligible resource oscillation. Using the
% resource avoids false cycle detections when grazer biomass is nearly zero.
cycle_relative_resource_amplitude = 1e-3;
cycle_absolute_resource_amplitude = 1e-4;

% Require the amplitudes in two consecutive late-time windows to agree.
settled_amplitude_ratio_bounds = [0.90,1.10];

%% Loop for filling the heatmaps
total_runs = nK*nW;
run_number = 0;

for w_index = 1:nW
    w = w_range(w_index);

    for K_index = 1:nK
        K = K_range(K_index);
        run_number = run_number+1;

        fprintf(['Run %d/%d: K = %.3f, w = %.3f\n'], ...
            run_number,total_runs,K,w);

        historyFcn = @(phi) history( ...
            phi,params,K,w,x_dim,J_dim,A_dim);

        try
            sol = dde23( ...
                @(phi,y,Z) dde(phi,y,Z,params,K,w), ...
                lag,historyFcn,phi_span,options);

            metrics = analyze_solution( ...
                sol,K,b,c_j,c_a, ...
                late_window_days,sample_step_days, ...
                extinction_threshold, ...
                cycle_relative_resource_amplitude, ...
                cycle_absolute_resource_amplitude, ...
                settled_amplitude_ratio_bounds);

            outcome_code(w_index,K_index) = metrics.outcome_code;

            if metrics.outcome_code == 1
                % Extinction is stored as NaN and drawn in white.
                grazer_amplitude(w_index,K_index) = NaN;
                cycle_period_days(w_index,K_index) = NaN;
                mean_stage_duration_days(w_index,K_index) = NaN;

            elseif metrics.outcome_code == 2
                % Persistent equilibrium: zero cycle amplitude and period.
                grazer_amplitude(w_index,K_index) = 0;
                cycle_period_days(w_index,K_index) = 0;
                mean_stage_duration_days(w_index,K_index) = ...
                    metrics.mean_stage_duration_days;

            elseif metrics.outcome_code == 3
                grazer_amplitude(w_index,K_index) = ...
                    metrics.grazer_amplitude;
                cycle_period_days(w_index,K_index) = ...
                    metrics.period_days;
                mean_stage_duration_days(w_index,K_index) = ...
                    metrics.mean_stage_duration_days;
            end

            fprintf(['  %s; amplitude = %.4g; period = %.4g days; ', ...
                'mean stage duration = %.4g days\n'], ...
                metrics.outcome_label, ...
                grazer_amplitude(w_index,K_index), ...
                cycle_period_days(w_index,K_index), ...
                mean_stage_duration_days(w_index,K_index));

        catch solver_error
            warning('K=%.3f, w=%.3f failed: %s', ...
                K,w,solver_error.message);
            outcome_code(w_index,K_index) = 0;
        end
    end

    % Save a checkpoint after every completed row of the heatmap.
    save('combined_model_heatmap_checkpoint.mat', ...
        'K_range','w_range','grazer_amplitude', ...
        'cycle_period_days','mean_stage_duration_days', ...
        'outcome_code','w_index');
end

%% Save completed results
save('combined_model_heatmap_results.mat', ...
    'K_range','w_range','grazer_amplitude', ...
    'cycle_period_days','mean_stage_duration_days', ...
    'outcome_code');

%% Local functions

function output = dde(~,y,Z,params,K,w)
    % Unpack parameters
    [kappa_j,kappa_a,r_j,alpha,beta, ...
        chi_j,chi_a,a_h,P,q,c_j,~,e_j,b] = ...
        deal(params{:});

    % Construct final nondimensional parameters
    lambda = a_h/K;
    Pi = P/(q*K);
    rho = (e_j*c_j)/(w*b);

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

    if ~isfinite(Q) || Q <= 0 || ~isfinite(D) || D <= 0
        error('The solution left the feasible region.');
    end

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

function yhist = history(phi,params,K,w,x_dim,J_dim,A_dim)
    % Unpack parameters needed for the history
    [kappa_j,kappa_a,r_j,~,~,~,~, ...
        a_h,P,q,c_j,c_a,e_j,b] = deal(params{:});

    lambda = a_h/K;
    Pi = P/(q*K);
    rho = (e_j*c_j)/(w*b);

    % Construct constant biological histories
    x0 = x_dim/K;
    J0 = (J_dim*c_j)/(b*K);
    A0 = (A_dim*c_a)/(b*K);

    % Evaluate the constant-history development rate
    Q0 = (Pi-kappa_j*J0-kappa_a*A0)/x0;
    h_j0 = tanh(r_j*Q0);
    F0 = x0/(lambda+x0);
    D0 = h_j0*F0;

    % Q0>1 is the nondimensional form of dimensional Q0>q.
    if x0 <= 0 || J0 < 0 || A0 < 0 || Q0 <= 1 || D0 <= 0
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

function metrics = analyze_solution( ...
        sol,K,b,c_j,c_a,late_window_days,sample_step_days, ...
        extinction_threshold,relative_cycle_threshold, ...
        absolute_cycle_threshold,settled_ratio_bounds)

    %% Convert to chronological time and dimensional variables
    time_days = sol.y(4,:)/b;
    resource = K*sol.y(1,:);
    juvenile = (b*K/c_j)*sol.y(2,:);
    adult = (b*K/c_a)*sol.y(3,:);
    grazer = juvenile+adult;

    phi = sol.x;

    %% Remove duplicates and invalid points
    [time_days,unique_indices] = unique(time_days,'stable');
    resource = resource(unique_indices);
    grazer = grazer(unique_indices);
    phi = phi(unique_indices);

    finite_values = ...
        isfinite(time_days) ...
        & isfinite(resource) ...
        & isfinite(grazer) ...
        & isfinite(phi);

    time_days = time_days(finite_values);
    resource = resource(finite_values);
    grazer = grazer(finite_values);
    phi = phi(finite_values);

    if numel(time_days) < 20 || any(diff(time_days) <= 0)
        error('Insufficient chronological-time output.');
    end

    %% Uniform chronological-time sampling
    analysis_start = max( ...
        time_days(end)-2*late_window_days, ...
        time_days(1));

    uniform_time = analysis_start:sample_step_days:time_days(end);

    if numel(uniform_time) < 40
        error('The simulation is too short for late-time analysis.');
    end

    uniform_resource = interp1( ...
        time_days,resource,uniform_time,'pchip');
    uniform_grazer = interp1( ...
        time_days,grazer,uniform_time,'pchip');
    uniform_phi = interp1( ...
        time_days,phi,uniform_time,'pchip');

    final_start = uniform_time(end)-late_window_days;
    preceding_start = final_start-late_window_days;

    preceding = ...
        uniform_time >= preceding_start ...
        & uniform_time < final_start;

    final = uniform_time >= final_start;

    if nnz(preceding) < 20 || nnz(final) < 20
        error('Two complete late-time windows are unavailable.');
    end

    %% Extinction and amplitude metrics
    final_mean_grazer = mean(uniform_grazer(final));
    final_max_grazer = max(uniform_grazer(final));

    final_resource_amplitude = ...
        max(uniform_resource(final))-min(uniform_resource(final));

    preceding_resource_amplitude = ...
        max(uniform_resource(preceding))-min(uniform_resource(preceding));

    final_grazer_amplitude = ...
        max(uniform_grazer(final))-min(uniform_grazer(final));

    relative_resource_amplitude = ...
        final_resource_amplitude ...
        /max(abs(mean(uniform_resource(final))),1e-12);

    amplitude_ratio = ...
        final_resource_amplitude ...
        /max(preceding_resource_amplitude,1e-12);

    is_extinct = ...
        final_mean_grazer < extinction_threshold ...
        && final_max_grazer < 10*extinction_threshold;

    is_cycle = ...
        ~is_extinct ...
        && final_resource_amplitude >= absolute_cycle_threshold ...
        && relative_resource_amplitude >= relative_cycle_threshold ...
        && amplitude_ratio >= settled_ratio_bounds(1) ...
        && amplitude_ratio <= settled_ratio_bounds(2);

    %% Period estimate
    period_days = NaN;

    if is_cycle
        period_days = estimate_cycle_period( ...
            uniform_time(final),uniform_resource(final));

        if ~isfinite(period_days)
            is_cycle = false;
        end
    end

    %% Mean stage duration
    %
    % For the transformed system,
    % tau(t(phi)) = [T(phi)-T(phi-1)]/b.
    final_phi = uniform_phi(final);
    valid_delay = final_phi >= 1;

    if any(valid_delay) && ~is_extinct
        delayed_clock = deval( ...
            sol,final_phi(valid_delay)-1,4)/b;

        current_time = uniform_time(final);
        current_time = current_time(valid_delay);

        stage_duration = current_time-delayed_clock;
        stage_duration = stage_duration( ...
            isfinite(stage_duration) & stage_duration >= 0);

        if isempty(stage_duration)
            mean_stage_duration = NaN;
        else
            mean_stage_duration = mean(stage_duration);
        end
    else
        mean_stage_duration = NaN;
    end

    %% Store classification
    if is_extinct
        code = 1;
        label = 'extinction';
    elseif is_cycle
        code = 3;
        label = 'limit cycle';
    else
        code = 2;
        label = 'persistent non-cycle';
    end

    metrics.outcome_code = code;
    metrics.outcome_label = label;
    metrics.grazer_amplitude = final_grazer_amplitude;
    metrics.period_days = period_days;
    metrics.mean_stage_duration_days = mean_stage_duration;
end

function period = estimate_cycle_period(time,signal)
    period = NaN;

    if numel(signal) < 5
        return
    end

    % Identify local maxima without requiring findpeaks.
    peak_indices = find( ...
        signal(2:end-1) > signal(1:end-2) ...
        & signal(2:end-1) >= signal(3:end))+1;

    signal_amplitude = max(signal)-min(signal);

    if signal_amplitude <= 0
        return
    end

    % Retain peaks in the upper quarter of the oscillation.
    threshold = min(signal)+0.75*signal_amplitude;
    peak_indices = peak_indices( ...
        signal(peak_indices) >= threshold);

    if numel(peak_indices) >= 3
        period = median(diff(time(peak_indices)));
    end
end

function limits = finite_color_limits(data)
    finite_data = data(isfinite(data));

    if isempty(finite_data)
        limits = [0,1];
        return
    end

    lower_limit = min(finite_data);
    upper_limit = max(finite_data);

    % Keep zero mapped to the dark-blue first color.
    lower_limit = min(0,lower_limit);

    if upper_limit <= lower_limit
        upper_limit = lower_limit+1;
    end

    limits = [lower_limit,upper_limit];
end
