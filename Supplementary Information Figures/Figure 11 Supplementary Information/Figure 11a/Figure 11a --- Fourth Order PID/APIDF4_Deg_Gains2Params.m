function [Parameters, Feasible] = APIDF4_Deg_Gains2Params(Gains, Parameters, r, SupportingInput)
%% Extract Gains
Feasible = true;
K_P = Gains.K_P;
K_I = Gains.K_I;
K_D = Gains.K_D;
omega_0 = Gains.omega_0;

%% Compute Supporting Input
[u_bar, x] = SupportingInput(Parameters, r);
x_bar_1 = x(1);

%% Extract Controller Parameters
mu = Parameters.mu;
kappa = Parameters.kappa;
n = Parameters.n;

%% Compute delta
delta = K_P * (1 + (kappa/x_bar_1)^n);

%% Compute k
k = K_I;

%% Compute eta_1
alpha_1 = omega_0;

%% Compute mu_1
alpha_2 = omega_0*K_D;

%% Compute theta
theta = mu/r;

%% Store Computed Parameters
Parameters.theta = theta;
Parameters.k = k;
Parameters.delta = delta;
Parameters.alpha_1 = alpha_1;
Parameters.alpha_2 = alpha_2;
end

