%% Basic Bifurcation Analysis of Combined Model
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

%% Delays
lag = 1;

%% History
x_dim = 0.3125;
J_dim = 0.09;
A_dim = 0.01;

%% K Range
K_range = 0.1:0.05:3;

%% Phi tspans
phi_span = [0, 5000];

%% options
options = ddeset('RelTol', 1e-6, 'AbsTol', 1e-9);

%% Min-max storage

x_min_storage = nan(size(K_range));
x_max_storage = nan(size(K_range));

grazer_min_storage = nan(size(K_range));
grazer_max_storage = nan(size(K_range));

%% Bifurcation loop

for K_index = 1:numel(K_range)

    K = K_range(K_index);

    fprintf('Current K value: %.3f\n',K);

    historyFcn = @(phi) history( ...
        phi,params,K,x_dim,J_dim,A_dim);

    sol = dde23( ...
        @(phi,y,Z) dde(phi,y,Z,params,K), ...
        lag,historyFcn,phi_span,options);

    %% Convert to dimensional variables

    t_days = sol.y(4,:)/b;

    x_dimensional = K*sol.y(1,:);
    J_dimensional = (b*K/c_j)*sol.y(2,:);
    A_dimensional = (b*K/c_a)*sol.y(3,:);

    grazer_dimensional = ...
        J_dimensional+A_dimensional;

    %% Retain only the final 5000 chronological days

    late_window_days = 5000;

    late = t_days >= ...
        max(t_days(end)-late_window_days,t_days(1));

    if nnz(late) < 10
        warning( ...
            'Too few late-time points at K = %.3f.',K);
        continue
    end

    %% Store late-time extrema

    x_min_storage(K_index) = ...
        min(x_dimensional(late));

    x_max_storage(K_index) = ...
        max(x_dimensional(late));

    grazer_min_storage(K_index) = ...
        min(grazer_dimensional(late));

    grazer_max_storage(K_index) = ...
        max(grazer_dimensional(late));
end

figure;
tiledlayout(1,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

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
ylabel('$J+A\;(\mathrm{mg\,C\,L^{-1})}$');
title('Total grazer bifurcation diagram');
grid on
box on

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