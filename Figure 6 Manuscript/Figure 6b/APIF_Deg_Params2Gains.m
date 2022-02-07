function Gains = APIF_Deg_Params2Gains(Parameters, SupportingInput)
% u = k * z_1   - delta * x_L * (x_1/kappa)^n / (1 + (x_1/kappa)^n)

%% Extract Controller Parameters
mu = Parameters.mu;
theta = Parameters.theta;
k = Parameters.k;
eta = Parameters.eta;
kappa = Parameters.kappa;
delta = Parameters.delta;
n = Parameters.n;
r = mu / theta;

%% Fixed Point
X_bar = APIF_Deg_FixedPoint(Parameters, SupportingInput);
Z_bar_1 = X_bar(end-1);
Z_bar_2 = X_bar(end);
X_bar_1 = X_bar(1);

%% Compute Derivatives
sigma_1 = k;
sigma_3 = delta * n * (r/X_bar_1) * (X_bar_1/kappa)^n / (1 + (X_bar_1/kappa)^n)^2;
sigma_4 = delta * (X_bar_1/kappa)^n / (1 + (X_bar_1/kappa)^n);

%% Compute PID Gains
K_P = sigma_4;
K_I = sigma_1 * Z_bar_1 / (Z_bar_1 + Z_bar_2);
omega_c = eta * (Z_bar_1 + Z_bar_2);
K_F = sigma_1 / (eta * (Z_bar_1 + Z_bar_2));
K_S = theta;

%% Construct Gains Structure
Gains.K_P = K_P;
Gains.K_I = K_I;
Gains.omega_c = omega_c;
Gains.K_F = K_F;
Gains.K_S = K_S;

end

