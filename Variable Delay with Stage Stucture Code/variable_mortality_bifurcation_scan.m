%% Basic Bifurcation-Style Scan of the Variable-Mortality Combined Model
clear;
close all;
clc;

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Dimensional parameter values
P = 0.025;       % Total phosphorus in system (mg P L^-1)
b = 1.2;         % Maximum producer growth rate (d^-1)
q = 0.0038;      % Minimum producer P:C ratio
e_j = 0.5;       % Maximum juvenile production efficiency
e_a = 0.8;       % Maximum adult production efficiency
theta_j = 0.025; % Juvenile P:C ratio
theta_a = 0.03;  % Adult P:C ratio
a_h = 0.25;      % Half-saturation constant (mg C L^-1)
c_j = 0.5;       % Maximum juvenile ingestion rate (d^-1)
c_a = 0.81;      % Maximum adult ingestion rate (d^-1)

%% Developmental threshold and variable-mortality scalers
w = 3;
mu_j = 7.45e-3;
mu_a = 5.92e-6;

%% Nondimensional parameters
kappa_j = (theta_j*b)/(c_j*q);
kappa_a = (theta_a*b)/(c_a*q);
r_j = q/theta_j;
r_a = q/theta_a;
alpha = (c_j*e_a)/b;
beta = (c_a*e_a)/b;
chi_j = mu_j/(e_j*c_j*b);
chi_a = mu_a/(e_a*c_a*b);
rho = (e_j*c_j)/(w*b);

params = {kappa_j,kappa_a,r_j,r_a,alpha,beta,chi_j,chi_a, ...
    rho,a_h,P,q,c_j,c_a,b};

%% Delay and constant dimensional history
lag = 1;

x_dim = 0.9;
J_dim = 0.04;
A_dim = 0.04;

%% Scan settings
K_range = 0.1:0.05:3;
phi_span = [0,5000];
late_window_days = 5000;

options = ddeset( ...
    'RelTol',1e-6, ...
    'AbsTol',1e-9, ...
    'MaxStep',1);

%% Storage
x_min_storage = nan(size(K_range));
x_max_storage = nan(size(K_range));
grazer_min_storage = nan(size(K_range));
grazer_max_storage = nan(size(K_range));
completed = false(size(K_range));

%% Bifurcation-style parameter scan
for K_index = 1:numel(K_range)

    K = K_range(K_index);
    fprintf('Current K value: %.3f (%d of %d)\n', ...
        K,K_index,numel(K_range));

    try
        historyFcn = @(phi) variable_mortality_history( ...
            phi,params,K,x_dim,J_dim,A_dim);

        sol = dde23( ...
            @(phi,y,Z) variable_mortality_dde(y,Z,params,K), ...
            lag,historyFcn,phi_span,options);

        %% Convert to dimensional variables
        t_days = sol.y(4,:)/b;
        x_dimensional = K*sol.y(1,:);
        J_dimensional = (b*K/c_j)*sol.y(2,:);
        A_dimensional = (b*K/c_a)*sol.y(3,:);
        grazer_dimensional = J_dimensional+A_dimensional;

        %% Retain the final chronological-time window
        late_start = max(t_days(end)-late_window_days,t_days(1));
        late = t_days >= late_start;

        if nnz(late) < 10
            warning('Too few late-time points at K = %.3f.',K);
            continue
        end

        %% Store late-time extrema
        x_min_storage(K_index) = min(x_dimensional(late));
        x_max_storage(K_index) = max(x_dimensional(late));

        grazer_min_storage(K_index) = min(grazer_dimensional(late));
        grazer_max_storage(K_index) = max(grazer_dimensional(late));

        completed(K_index) = true;

    catch ME
        warning('K = %.3f was skipped: %s',K,ME.message);
    end
end

%% Plot using the same two-panel structure as the constant-death scan
figure('Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile
plot(K_range,x_max_storage,'r.','MarkerSize',10);
hold on
plot(K_range,x_min_storage,'r.','MarkerSize',10);

xlabel('$K\;(\mathrm{mg\,C\,L^{-1}})$');
ylabel('$\mathrm{Resource\ density\;(mg\,C\,L^{-1})}$');
title('Resource bifurcation diagram');
grid on
box on

nexttile
plot(K_range,grazer_max_storage,'k.','MarkerSize',10);
hold on
plot(K_range,grazer_min_storage,'k.','MarkerSize',10);

xlabel('$K\;(\mathrm{mg\,C\,L^{-1}})$');
ylabel('$J+A\;(\mathrm{mg\,C\,L^{-1}})$');
title('Total grazer bifurcation diagram');
grid on
box on

tilte("$w = 3, \mu_j = ")

title_text = ...
    "Variable mortality: $w=" + compose("%.2f",w) + ...
    "$, $\mu_j=" + compose("%.3g",mu_j) + ...
    "$, $\mu_a=" + compose("%.3g",mu_a) + "$";

sgtitle(title_text, 'Interpreter', 'latex');

save('variable_mortality_bifurcation_scan_results.mat', ...
    'K_range','w','mu_j','mu_a', ...
    'x_min_storage','x_max_storage', ...
    'grazer_min_storage','grazer_max_storage','completed');

%% DDE system
function output = variable_mortality_dde(y,Z,params,K)

    [kappa_j,kappa_a,r_j,r_a,alpha,beta,chi_j,chi_a, ...
        rho,a_h,P,q,~,~,~] = deal(params{:});

    % State order: [x; J; A; T; R]
    output = zeros(5,1);

    x = y(1);
    J = y(2);
    A = y(3);

    lambda = a_h/K;
    Pi = P/(q*K);

    if x <= 0
        error('Resource biomass became nonpositive.');
    end

    Q = (Pi-kappa_j*J-kappa_a*A)/x;

    if Q <= 0
        error('The resource P:C ratio became nonpositive.');
    end

    h_j = tanh(r_j*Q);
    h_a = tanh(r_a*Q);
    F = x/(lambda+x);
    D_j = h_j*F;
    D_a = h_a*F;

    if D_j <= 0 || D_a <= 0
        error('A developmental or production rate became nonpositive.');
    end

    %% Through-stage survival from the cumulative integral state
    survival_integral = y(5)-Z(5,1);

    integral_tolerance = 1e-10*max( ...
        [1,abs(y(5)),abs(Z(5,1))]);

    if survival_integral < -integral_tolerance
        error('The survival integral became negative: %.6e.', ...
            survival_integral);
    end

    survival_integral = max(survival_integral,0);
    S_phi = exp(-(chi_j/rho)*survival_integral);

    %% Population equations
    output(1) = ...
        (x/(rho*D_j))*(1-x)*(1-1/Q) ...
        -(J+A)/(rho*h_j);

    output(2) = ...
        (alpha/rho)*(A-Z(3,1)*S_phi) ...
        -(chi_j/(rho*D_j^2))*J;

    output(3) = ...
        (beta/rho)*Z(3,1)*S_phi ...
        -(chi_a/(rho*D_j*D_a))*A;

    %% Chronological clock and survival-integral accumulator
    output(4) = 1/(rho*D_j);
    output(5) = 1/D_j^2;
end

%% Constant history
function yhist = variable_mortality_history( ...
        phi,params,K,x_dim,J_dim,A_dim)

    [kappa_j,kappa_a,r_j,~,~,~,~,~,rho,a_h,P,q,c_j,c_a,b] = ...
        deal(params{:});

    lambda = a_h/K;
    Pi = P/(q*K);

    x0 = x_dim/K;
    J0 = (J_dim*c_j)/(b*K);
    A0 = (A_dim*c_a)/(b*K);

    Q0 = (Pi-kappa_j*J0-kappa_a*A0)/x0;
    h_j0 = tanh(r_j*Q0);
    F0 = x0/(lambda+x0);
    D_j0 = h_j0*F0;

    if x0 <= 0 || Q0 <= 0 || D_j0 <= 0
        error('The selected history is not biologically feasible.');
    end

    nPhi = numel(phi);

    % Linear T and R histories make their delayed differences consistent
    % with the constant biological history on phi in [-1,0].
    yhist = [ ...
        x0*ones(1,nPhi);
        J0*ones(1,nPhi);
        A0*ones(1,nPhi);
        phi/(rho*D_j0);
        phi/(D_j0^2)
    ];
end
