function [Parameters, Feasible] = APIF_Deg_Gains2Params(Gains, Parameters, r, SupportingInput)
%% Extract Gains
Feasible = true;
K_P = Gains.K_P;
K_I = Gains.K_I;
omega_c = Gains.omega_c;

%% Compute Supporting Input
[u_bar, x] = SupportingInput(Parameters, r);
x_bar_1 = x(1);

%% Extract Controller Parameters
mu = Parameters.mu;
kappa = Parameters.kappa;
n = Parameters.n;

%% Compute delta
delta = K_P * (1 + (kappa/x_bar_1)^n);

%% Compute eta
eta = K_I*omega_c / (u_bar + r*K_P);

%% Compute k
k = (omega_c/mu/2) * (u_bar + r*K_P) * (1 - sqrt(1 - 4*mu*K_I/(u_bar + K_P*r)/omega_c ));
if imag(k)~=0
    k = NaN;
    Feasible = false;
end

%% Compute theta
theta = mu/r;

%% Store Computed Parameters
Parameters.delta = delta;
Parameters.eta = eta;
Parameters.k = k;
Parameters.theta = theta;
end

