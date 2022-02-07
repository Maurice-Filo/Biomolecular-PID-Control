function [Parameters, Feasible] = APIDF4H_Deg_Gains2Params(Gains, Parameters, r, SupportingInput)
% u = k * z_1 + k_0 * x_L - (alpha_1*z_3 + alpha_2*x_L) * (x_1/kappa)^n / (1 + (x_1/kappa)^n)

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

%% Compute k_0
k_0 = abs(K_P);

%% Compute k
k = K_I;

%% Compute alpha_1
alpha_1 = omega_0;

%% Compute alpha_2
alpha_2 = omega_0*K_D;

%% Compute theta
theta = mu/r;

%% Store Computed Parameters
Parameters.theta = theta;
Parameters.k = k;
Parameters.k_0 = k_0;
Parameters.alpha_1 = alpha_1;
Parameters.alpha_2 = alpha_2;

end

