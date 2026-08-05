%% Resource Dependent Death Combined Model Simulations
% This script will simulate our resource dependent model for various values
% of our carrying capacity K outputing solutions, Q plots, and stage
% duration vs period graphs. 

clear;
close all;
clc;

%%  Parameter Values (Taken from Tomas)
P = 0.025;       %Total phosphorous in system
b = 1.2;         %Producer max growth rate
q = 0.0038;      %Producer min P:C
e_j = 0.5;       %Juevenile max production efficiency
e_a = 0.8;       %Adult max production efficiency
theta_j = 0.025; %Juvenile constant P:C ratio
theta_a = 0.03;  %Adult constant P:C ratio
a_h = 0.25;      %Juvenile and adult half-saturation constant
c_j = 0.5;       %Juvenile max ingestion rate
c_a = 0.81;      %Adult max ingestion rate 

%% Stage Duration Parameter
w = 1.7;

%% Mortality Scalers
mu_j = 0.00128;
mu_a = 0.000560;

%% Low-Pass Filter Approximation Terms
nu = 1;
n = 5;
zeta = 1000;

filter_params = {nu, n, zeta};

%% Non-Dimensional Parameters
kappa_j = (theta_j*b)/(c_j*q);
kappa_a = (theta_a*b)/(c_a*q);
r_j = q/theta_j;
r_a = q/theta_a;
alpha = (c_j*e_a)/b;
beta = (c_a*e_a)/b;
chi_j = mu_j/(e_j*c_j*b);
chi_a = mu_a/(e_a*c_a*b);
rho = (e_j*c_j)/(w*b);

params = {kappa_j, kappa_a, r_j, r_a, alpha, beta, chi_j, chi_a, rho, a_h, P, q, c_j, c_a, b};

%% K values
K_1 = 1.5;

%% Constant dimensional history
x_dim = 1.08;
J_dim = 0.04;
A_dim = 0.04;

%% Constant non-dimensional history
x0 = x_dim/K_1;
J0 = (J_dim*c_j)/(b*K_1);
A0 = (A_dim*c_a)/(b*K_1);

%% Helper functions at the initial history

lambda_0 = a_h/K_1;
Pi_0 = P/(q*K_1);

Q0 = ...
    (Pi_0-kappa_j*J0-kappa_a*A0)/x0;

h_j0 = tanh(r_j*Q0);
F0 = x0/(lambda_0+x0);

D_j0 = h_j0*F0;

if x0 <= 0 || Q0 <= 0 || D_j0 <= 0
    error('The selected initial history is not feasible.');
end

%% History function

% State order: [x; J; A; t; R]
historyFcn_1 = @(phi) [
    x0*ones(1,numel(phi));
    J0*ones(1,numel(phi));
    A0*ones(1,numel(phi));
    phi/(rho*D_j0);
    phi/(D_j0^2)
];

%% Delay and developmental-time span
lag = 1;
phi_span = [0, 500];

%% options
options = ddeset( ...
    'RelTol',1e-8, ...
    'AbsTol',1e-11, ...
    'MaxStep',0.01);

%% solutions
sol_K_1 = dde23(@(phi,y,Z) dde(phi, y, Z, params, K_1), lag, historyFcn_1, phi_span, options);

%% Convert solution to dimensional variables
phi_1 = sol_K_1.x;

x_1_non_dim = sol_K_1.y(1,:);
J_1_non_dim = sol_K_1.y(2,:);
A_1_non_dim = sol_K_1.y(3,:);
T_1_non_dim = sol_K_1.y(4,:);

% Convert chronological time to days
t_1 = T_1_non_dim/b;

% Convert biomass variables to dimensional units
x_1 = K_1*x_1_non_dim;
J_1 = (b*K_1/c_j)*J_1_non_dim;
A_1 = (b*K_1/c_a)*A_1_non_dim;

grazer_1 = J_1+A_1;

% Avoid division by zero when calculating adult proportion
adult_proportion_1 = A_1./max(grazer_1,eps);

%% Plot population solution

figure('Color','w','Name','Population time series');
clf

yyaxis left

plot(t_1,grazer_1, ...
    'k-','LineWidth',1.5);
hold on

plot(t_1,x_1, ...
    'k--','LineWidth',1.5);

ylabel('$\mathrm{Density\;(mg\,C\,L^{-1})}$');
xlim([0, 500]);

yyaxis right

plot(t_1,adult_proportion_1, ...
    'r-','LineWidth',1.5);

ylabel('$\mathrm{Adult\ proportion}$');
ylim([0 1]);

ax = gca;
ax.YAxis(2).Color = 'r';

xlabel('$t\;(\mathrm{days})$');

title(sprintf('$K=%.2f$, $w=%.2f$',K_1,w));

legend({ ...
    '$\mathrm{Grazer}\;(J+A)$', ...
    '$\mathrm{Producer}$', ...
    '$\mathrm{Adult\ proportion}$'}, ...
    'Location','best');

box on
grid on

%% Stage duration as a function of chronological time

% t(phi-1) is available from the solution for phi >= 1
valid_delay = phi_1 >= 1;

phi_delay = phi_1(valid_delay);
T_now = T_1_non_dim(valid_delay);

% Evaluate T(phi-1)
T_delayed = deval( ...
    sol_K_1,phi_delay-1,4);

% Convert to chronological days
t_delay_days = T_now/b;

tau_1_days = ...
    (T_now-T_delayed)/b;

%% Cycle period as a function of time

% Ignore extremely small oscillations
grazer_range = range(grazer_1);

minimum_prominence = max( ...
    0.01*grazer_range,1e-6);

% Detect grazer maxima in chronological time
[~,peak_times] = findpeaks( ...
    grazer_1,t_1, ...
    'MinPeakProminence',minimum_prominence);

if numel(peak_times) >= 2

    % Period of every completed cycle
    cycle_period = diff(peak_times);

    % Associate each period with the midpoint between its two peaks
    period_time = ...
        0.5*(peak_times(1:end-1)+peak_times(2:end));

else

    cycle_period = [];
    period_time = [];

    warning('Not enough peaks were detected to calculate cycle periods.');

end

%% Stage-duration and time-dependent period plot

figure('Color','w','Name','Stage timing');
clf

plot(t_delay_days,tau_1_days, ...
    'k-','LineWidth',1.5);
hold on

if ~isempty(cycle_period)
    plot(period_time,cycle_period, ...
        '--', ...
        'Color',[0.65 0.65 0.65], ...
        'LineWidth',1.5);
end

xlabel('$t\;(\mathrm{days})$');
xlim([0, 500]);

ylabel( ...
    '$\mathrm{Stage\ duration\ and\ cycle\ period\;(days)}$');

title(sprintf('$K=%.2f$, $w=%.2f$',K_1,w));

if ~isempty(cycle_period)
    legend({ ...
        '$\tau(t)$', ...
        '$\mathrm{Cycle\ period}$'}, ...
        'Location','best');
else
    legend('$\tau(t)$','Location','best');
end

plot_values = tau_1_days(:);

if ~isempty(cycle_period)
    plot_values = [
        plot_values;
        cycle_period(:)
    ];
end

upper_limit = max(plot_values,[],'omitnan');

if isfinite(upper_limit) && upper_limit > 0
    ylim([0,1.1*upper_limit]);
end

box on

function output = dde(~,y,Z,params,K)

    %% Unpack parameters

    [kappa_j, kappa_a, r_j, r_a, alpha, beta, ...
        chi_j, chi_a, rho, a_h, P, q, ~, ~, ~] = ...
        deal(params{:});

    % Five states: [x; J; A; t; R]
    output = zeros(5,1);

    %% Extract current states

    x = y(1);
    J = y(2);
    A = y(3);

    %% Define dimensionless parameters and helper functions

    lambda = a_h/K;
    Pi = P/(q*K);

    Q = ...
        (Pi-kappa_j*J-kappa_a*A)/x;

    h_j = tanh(r_j*Q);
    h_a = tanh(r_a*Q);

    F = x/(lambda+x);

    D_j = h_j*F;
    D_a = h_a*F;

    %% Check biological feasibility

    if x <= 0
        error('Resource biomass became nonpositive.');
    end

    if Q <= 0
        error('The resource P:C ratio became nonpositive.');
    end

    if D_j <= 0 || D_a <= 0
        error('A developmental or production rate became nonpositive.');
    end

    %% Calculate survival directly from the accumulator state

    % R(phi)-R(phi-1) equals the survival integral
    survival_integral = ...
        y(5)-Z(5,1);

    % Scale tolerance to the size of the cumulative states
    integral_tolerance = ...
        1e-10*max([ ...
        1,abs(y(5)),abs(Z(5,1))]);

    if survival_integral < -integral_tolerance
        error( ...
            ['The calculated survival integral became ', ...
             'negative: %.6e'], ...
            survival_integral);
    end

    % Remove only tiny negative roundoff errors
    survival_integral = ...
        max(survival_integral,0);

    S_phi = exp(-(chi_j/rho)*survival_integral);

    %% Population equations

    output(1) = (x/(rho*D_j))*(1-x) *(1-1/Q) -(J+A)/(rho*h_j);

    output(2) = (alpha/rho) *(A-Z(3,1)*S_phi)-(chi_j/(rho*D_j^2))*J;

    output(3) = (beta/rho) *Z(3,1)*S_phi -(chi_a/(rho*D_j*D_a))*A;

    %% Chronological-time equation

    output(4) = 1/(rho*D_j);

    %% Cumulative survival-integral equation

    output(5) =  1/D_j^2;
end
