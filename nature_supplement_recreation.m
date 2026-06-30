%% Recreating Figure 1 from Nature Supplement 

%Daphnia = Black
%Algea = Gray
%Cycle period = dashed line
%Stage duration = solid line

%%  Parameter Values
q = 1;
I_J = 5.23 * 10^(-3);
I_A = 1.91 * 10^(-2);
f_h = 0.164;
chi = 0.77;
gamma = 1.51*10^(-3);
sigma_J = 0.49;
sigma_A = 0.43;
mu_J = 1.09*10^(-5);
mu_A = 5.92*10^(-6);
w = 4.8 * 10^(-3);

%%  Non-dimensional Parameter Values

alpha_J = (chi*I_J*sigma_A)/(gamma*q);
alpha_A = (chi*I_A*sigma_A)/(gamma*q);
nu_J = (mu_J)/(q*I_J*sigma_J);
nu_A = (mu_A)/(q*I_A*sigma_A);
rho = (I_J*sigma_J)/(w*q);

params = {alpha_J, alpha_A, nu_J, nu_A, rho, f_h};

%% Low-pass filter params
nu = 1;
n = 5;
zeta = 1000;

filter_params = {nu, n, zeta};

%% delays
lag = 1;

%% history
Q_initial = 0.1*ones(n, 1);
history = [27/500; 0.1; 0.1; Q_initial; 0];

%% phi_spans
phi_span_1 = [0, 250];
phi_span_2 = [0, 1000];

%% k values
k_1 = 0.6;
k_2 = 0.3;

%% options
options = ddeset('RelTol', 1e-6, 'AbsTol', 1e-9);

%% calculating solutions
sol_k_1 = dde23(@(phi,y,Z) dde(phi,y,Z,params, filter_params, k_1), lag, history, phi_span_1, options);
sol_k_2 = dde23(@(phi,y,Z) dde(phi,y,Z,params, filter_params, k_2), lag, history, phi_span_2, options);

%% Grabbing outputs
f_1 = sol_k_1.y(1,:);
j_1 = sol_k_1.y(2,:);
a_1 = sol_k_1.y(3,:);
t_1 = sol_k_1.y(9,:);

%% Grabbing input phi and calculating tau(t) 
% Collect X-axis of dde which has units of Phi
phi_vals_1 = sol_k_1.x;

% Note the 'valid' index values of phi which have value greater than 1
idx_1 = phi_vals_1 >= 1;

% Collect the 'valid' phi and corresponding time values. This corresponds
% to grabbing phi and t(phi)
phi_use_1 = phi_vals_1(idx_1);
t_now_1   = t_1(idx_1);

%Now with this 'valid' phi vector, evauate at all phi - 1 to eventually
%grab the correspndong t(phi-1_
Y_delay_1 = deval(sol_k_1, phi_use_1 - 1);
t_past_1  = Y_delay_1(9,:);

%compute the delay
tau_1 = t_now_1 - t_past_1;

%% Dimensionalization 
F_1 = f_1 * k_1;
J_1 = j_1 * (q*k_1/I_J);
A_1  = a_1 * (q*k_1/I_A);

%% Finding period over time
D_1 = F_1 + J_1;
[pks_1, peak_times_1] = findpeaks(D_1, t_1, ...
    'MinPeakDistance', 5);
period_1 = diff(peak_times_1);
period_time_1 = peak_times_1(2:end);

%% Plotting for k = 0.3;
figure(1);
plot(t_1, J_1 + A_1, 'k', 'LineWidth', 1.5) 
hold on
plot(t_1, F_1*500, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5)
xlabel('Days');
ylabel('Daphnia (#/L) & Food Density (mgC/L x500)');
title('(a) with k = 0.3')
xlim([0 250]);
ylim([0 50]);

figure(2);
plot(t_now_1, tau_1, 'k', 'LineWidth', 1.5)
hold on
plot(period_time_1, period_1, 'Color', [0.6 0.6 0.6], 'LineStyle', '--')
xlabel('Days')
ylabel('\tau_1(t) Stage Duration & Cycle Period (days)')
title('(b) with k = 0.3')
xlim([0 250])

%% Grabbing outputs
f_2 = sol_k_2.y(1,:);
j_2 = sol_k_2.y(2,:);
a_2 = sol_k_2.y(3,:);
t_2 = sol_k_2.y(9,:);

%% Grabbing input variable phi and calucalting tau(t)
phi_value_2 = sol_k_2.x;

idx_2 = phi_value_2 >= 1;

phi_use_2 = phi_value_2(idx_2);
t_now_2 = t_2(idx_2);

Y_delay_2 = deval(sol_k_2, phi_use_2 - 1);
t_past_2 = Y_delay_2(9,:);

tau_2 = t_now_2 - t_past_2;

%% Dimensionalization 
F_2 = f_2 * k_2;
J_2 = j_2 * (q*k_2/I_J);
A_2  = a_2 * (q*k_2/I_A);

%% Finding period over time
D_2 = F_2 + J_2;
[pks_2, peak_times_2] = findpeaks(D_2, t_2, ...
    'MinPeakDistance', 5);
period_2 = diff(peak_times_2);
period_time_2 = peak_times_2(2:end);

%% Plotting for k = 0.6
figure(3);
plot(t_2, J_2 + A_2, 'k', 'LineWidth', 1.5 );
hold on
plot(t_2,F_2*500, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5 )
xlabel('Days');
ylabel('Daphnia (#/L) & Food Density (mgC/L x500)');
title('(c) with k = 0.6')
xlim([0 1000]);
ylim([0 50]);

figure(4)
plot(t_now_2, tau_2, 'k', 'LineWidth', 1.5)
hold on
plot(period_time_2, period_2, 'Color', [0.6 0.6 0.6], 'LineStyle', '--')
xlabel('Days')
ylabel('\tau_1(t) Stage Duration (days)')
title('(d) with k = 0.6')
xlim([0 1000])

%% delay function
function output = dde(phi,y,Z,params, filter_params, k)
    %unpack params
    [alpha_J, alpha_A, nu_J, nu_A, rho, f_h] = params{:};
    lambda = f_h/k;
    [nu, n, zeta] = filter_params{:};

    %preallocate derivative vector 
    output = zeros(size(y));

    %extract variables
    F = y(1); 
    J = y(2);
    A = y(3);

    % building dQ
    % Q_1 equation
    output(4) = -nu*y(4) + (1/zeta)*((F + lambda)/F)^2;
    
    % dQ_i equations for i > 1
    for i = 2:n
        output(3+i) = -nu*y(3+i) + y(3+i-1);
    end

    weighted_sum = 0;

    for i = 1:n
        weighted_sum = weighted_sum + nu^(i-1)*(y(3+i) - Z(3+i,1));
    end

    inside = weighted_sum + nu^n*y(3+n);

    S_phi = exp(-zeta*nu_J/rho * inside);

    %poulation eqs
    output(1) = (1/rho)*((1-y(1))*(y(1)+lambda)-y(2)-y(3));
    output(2) = (alpha_J/rho)*(y(3)-Z(3,1)*S_phi) - (nu_J/rho)*((F+lambda)/F)^2*J;
    output(3) = (alpha_A/rho)*Z(3,1)*S_phi-(nu_A/rho)*((F+lambda)/F)^2*A;

    %time eq
    time_index = 4+n;
    output(time_index) = (F + lambda)/(rho*F);
end

