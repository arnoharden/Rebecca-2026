% Recreating Figure S1 transient dynamics

% Parameters used for dimensional scaling and plotting
k    = 0.6;        % resource carrying capacity K, mgC/L
q    = 1.0;        % max per-capita resource growth rate, day^-1
IJ   = 5.23e-3;   % max juvenile ingestion rate
IA   = 1.91e-2;   % max adult ingestion rate

% a fixed delay of length 1 in phi-space.
omega = 1;

% Integration range in phi-space.
phispan = [0 250];

options = ddeset('RelTol', 1e-6, 'AbsTol', 1e-9);


% First is history for K=0.3, second is history for K=0.6
history = [0.5; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.0];
% history = [28/500; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.0];

%% Solver
sol = dde23(@(phi,y,Z) ddefun(phi,y,Z,k), omega, history, phispan, options);

% y(9) stores real time in days, computed from dt/dphi.
t_days = sol.y(9,:);
t_days = t_days - t_days(1);

% re-dimensionalize
F_dim = sol.y(1,:) * k;
J_dim = sol.y(2,:) * (q*k/IJ);
A_dim = sol.y(3,:) * (q*k/IA);

%% Plotting

figure; hold on;

% Scale algae density x500
plot(t_days, F_dim*500, '-', ...
    'Color', [0.75 0.75 0.75], ...
    'LineWidth', 1.2, ...
    'DisplayName', 'Food density x500');

% Total Daphnia
plot(t_days, J_dim + A_dim, '-', ...
    'Color', [0 0 0], ...
    'LineWidth', 1.2, ...
    'DisplayName', 'Total Daphnia');

xlabel('Days');
ylabel('Daphnia (#/L) & Food Density (mgC/L; x500)');
xlim([0 250]);
ylim([0 50]);

legend('Location','northeast','Box','off');
box on;


%% Stage duration vs cycle period graphs

% sol.x is the independent variable phi.
phi_sol = sol.x;

% sol.y(9,:) is real time t(phi), in days.
t_days = sol.y(9,:);

%% Stage Duration
tau_days = zeros(size(phi_sol));

for i = 1:length(phi_sol)
    phi_target = phi_sol(i) - 1;

    if phi_target >= phi_sol(1)
        % Interpolate real time at phi - 1.
        t_past = interp1(phi_sol, t_days, phi_target);

        % Stage duration in days.
        tau_days(i) = t_days(i) - t_past;
    else
        tau_days(i) = NaN;
    end
end

total_daphnia = J_dim + A_dim;

[peak_vals, peak_t] = findpeaks(total_daphnia, t_days, 'MinPeakProminence', 0.2);

% Cycle periods are differences between successive peak times.
periods = diff(peak_t);

% Plot each period at the time of the second peak in the pair.
period_t = peak_t(2:end);

% Plot stage duration and cycle period.
figure; hold on;

% Solid black stage duration
plot(t_days, tau_days, '-', ...
    'Color', [0 0 0], ...
    'LineWidth', 1.2, ...
    'DisplayName', 'Juvenile stage duration');

% Dashed light gray cycle period
plot(period_t, periods, '--', ...
    'Color', [0.65 0.65 0.65], ...
    'LineWidth', 1.2, ...
    'DisplayName', 'Cycle period');

xlabel('Days');
ylabel('Stage Duration & Cycle Period (days)');
ylim([0 40]);
xlim([0 250]);

legend('Location','southeast','Box','off');
box on;

%% Function

% Model equations in transformed phi-time.
function dydt = ddefun(phi, y, Z, k)
    A_delay  = Z(3,1);
    Q1_delay = Z(4,1);
    Q2_delay = Z(5,1);
    Q3_delay = Z(6,1);
    Q4_delay = Z(7,1);
    Q5_delay = Z(8,1);

    % Dimensional parameter values.
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

    lambda = fh / k;

    alphaJ = (chi / gamma) * (IJ * sigA / q);
    alphaA = (chi / gamma) * (IA * sigA / q);

    etaJ = muJ / (q * IJ * sigJ);
    etaA = muA / (q * IA * sigA);

    rho = (IJ * sigJ) / (w * q);

    % Low-pass filter approximation settings.
    % The supplement uses nu = 1, n = 5, zeta = 1000.
    nu   = 1;
    zeta = 1000;

    % Survival probability approximation
    S = exp(-zeta*(etaJ/rho) * (y(4) - Q1_delay + nu   * (y(5) - Q2_delay) + nu^2 * (y(6) - Q3_delay) + nu^3 * (y(7) - Q4_delay) + nu^4 * (y(8) - Q5_delay) + nu^5 * y(8)));

    dydt = [
        (1/rho)*((1-y(1))*(y(1)+lambda) - y(2) - y(3));
        (alphaJ/rho)*(y(3) - A_delay*S) - (etaJ/rho)*((y(1)+lambda)/y(1))^2 * y(2);
        (alphaA/rho)*A_delay*S - (etaA/rho)*((y(1)+lambda)/y(1))^2 * y(3);
        -nu*y(4) + (1/zeta)*((y(1)+lambda)/y(1))^2;
        -nu*y(5) + y(4);
        -nu*y(6) + y(5);
        -nu*y(7) + y(6);
        -nu*y(8) + y(7);
        (y(1)+lambda) / (rho*y(1));
    ];
end