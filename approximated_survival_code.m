%% File containing dde23 function and history for approximated survival integral model
% To be used for bifurcation analysis

%% Low-Pass Filter Approximation Terms
nu = 1;
n = 5;
zeta = 1000;

filter_params = {nu, n, zeta};

%% Constant history
x_dim = 1.08;
J_dim = 0.04;
A_dim = 0.04;

%% Constant biological history
x0 = x_dim/K_1;
J0 = (J_dim*c_j)/(b*K_1);
A0 = (A_dim*c_a)/(b*K_1);

%% Development rate associated with the constant history

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

%% Consistent constant filter histories

z0 = 1/(zeta*D_j0^2);

L_initial = zeros(n,1);

for i = 1:n
    L_initial(i) = z0/nu^i;
end

%% Complete constant history vector

history = [
    x0;
    J0;
    A0;
    0;
    L_initial
];


function output = dde(phi,y,Z,params, filter_params, K)

    %unpack params
    [kappa_j, kappa_a, r_j, r_a, alpha, beta, chi_j, ...
    chi_a, rho, a_h, P, q, ~, ~, ~] = deal(params{:});

    [nu, n, zeta] = deal(filter_params{:});

    %preallocate derivative vector 
    output = zeros(size(y));

    %extract variables
    x = y(1); 
    J = y(2);
    A = y(3);
    
    %define extra parameters and functions
    lambda = a_h/K;
    Pi = P/(q*K);
    Q = (Pi - kappa_j*J - kappa_a*A)/x;
    h_j = tanh(r_j * Q); 
    h_a = tanh(r_a * Q);
    D_j = h_j*(x/(lambda + x));
    D_a = h_a*(x/(lambda + x));

    %% Survival approximation
    % L_1
    output(5) = -nu*y(5) + (1/zeta)*(1/D_j^2);
    
    % L_2 through L_n
    for i = 2:n
        current_index = 4 + i;
        previous_index = current_index - 1;
    
        output(current_index) = ...
            -nu*y(current_index) + y(previous_index);
    end
    
    % Sum of delayed filter differences
    weighted_sum = 0;
    
    for i = 1:n
        filter_index = 4 + i;
    
        weighted_sum = weighted_sum ...
            + nu^(i-1)*( ...
            y(filter_index) - Z(filter_index,1));
    end
    
    final_filter_index = 4 + n;
    
    inside = weighted_sum ...
        + nu^n*y(final_filter_index);
    
    S_phi = exp(-(zeta*chi_j/rho)*inside);

    if ~isfinite(S_phi)
        error('Nonfinite survival probability at K = %.4f.',K);
    end
    
    if S_phi > 1 + 1e-8 || S_phi < 0
        warning( ...
            'Survival probability outside [0,1]: S = %.6g', ...
            S_phi);
    end
    
    %% population equations
    output(1) = (x/(rho*D_j))*(1-x)*(1-(1/Q)) - (J+A)/(rho*h_j);
    output(2) = (alpha/rho)*(A - Z(3,1)*S_phi) - (chi_j/(rho*D_j^2))*J;
    output(3) = (beta/rho)*Z(3,1)*S_phi - (chi_a/(rho*D_j*D_a))*A;
    
    %% time equation
    output(4) = 1/(rho * D_j);
end