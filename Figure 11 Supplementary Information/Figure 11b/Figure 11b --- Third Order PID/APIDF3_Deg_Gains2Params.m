function [Parameters, Feasible] = APIDF3_Deg_Gains2Params(Gains, Parameters, r, SupportingInput)
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
alpha_0 = Parameters.alpha_0;
kappa = Parameters.kappa;
kappa_0 = Parameters.kappa_0;
n = Parameters.n;
n_0 = Parameters.n_0;

%% Compute delta
delta = K_D*omega_0 * (1 + (kappa/x_bar_1)^n);

%% Compute k
k = K_I;

%% Compute gamma_0
gamma_0 = omega_0;

%% Compute delta_0
delta_0 = omega_0 * (K_D*omega_0 - K_P) * (r / n_0 / alpha_0) * (1 + (kappa/x_bar_1)^n) * (1 + (r/kappa_0)^n_0)^2 / (r/kappa_0)^n_0;
if delta_0 < 0
    Feasible = false;
end

%% Compute theta
theta = mu/r;

%% Store Computed Parameters
Parameters.theta = theta;
Parameters.k = k;
Parameters.delta = delta;
Parameters.delta_0 = delta_0;
Parameters.gamma_0 = gamma_0;
end

