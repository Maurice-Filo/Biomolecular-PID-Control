function [Parameters, Feasible] = APIDF2_Deg_Gains2Params(Gains, Parameters, r, SupportingInput)
%% Extract Gains
Feasible = true;
K_P = Gains.K_P;
K_I = Gains.K_I;
K_D = Gains.K_D;
omega_c = Gains.omega_c;

%% Compute Supporting Input
[u_bar, x] = SupportingInput(Parameters, r);
x_bar_1 = x(1);

%% Extract Controller Parameters
mu = Parameters.mu;
kappa = Parameters.kappa;
n = Parameters.n;

%% Compute delta
delta = K_D*omega_c * (1 + (kappa/x_bar_1)^n);

%% Compute eta
T = u_bar + r*K_D*omega_c;
eta = K_I*omega_c / T;

%% Compute k
k = (omega_c/mu/2) * (u_bar + r*K_P - sqrt((u_bar + r*K_P)^2 - 4*mu*K_I*(u_bar + K_D*omega_c*r)/omega_c ));
if imag(k)~=0
    k = NaN;
    Feasible = false;
end

%% Compute beta
beta = omega_c * (K_D*omega_c - K_P) / k;
if beta < 0
    beta = NaN;
    Feasible = false;
end

%% Compute theta
theta = mu/r + beta;

%% Store Computed Parameters
Parameters.delta = delta;
Parameters.eta = eta;
Parameters.k = k;
Parameters.beta = beta;
Parameters.theta = theta;
end

