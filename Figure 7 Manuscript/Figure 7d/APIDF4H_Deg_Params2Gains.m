function Gains = APIDF4H_Deg_Params2Gains(Parameters, SupportingInput)
% u = k * z_1 + k_0 * x_L - (alpha_1*z_3 + alpha_2*x_L) * (x_1/kappa)^n / (1 + (x_1/kappa)^n)

%% Extract Controller Parameters
mu = Parameters.mu;
theta = Parameters.theta;
k = Parameters.k;
k_0 = Parameters.k_0;
mu_0 = Parameters.mu_0;
alpha_1 = Parameters.alpha_1;
alpha_2 = Parameters.alpha_2;
eta = Parameters.eta;
eta_0 = Parameters.eta_0;
kappa = Parameters.kappa;
n = Parameters.n;
delta_c = Parameters.delta_c;
r = mu / theta;

%% Fixed Point
X_bar = APIDF4H_Deg_FixedPoint(Parameters, SupportingInput);
Z_bar_1 = X_bar(end-3);
Z_bar_2 = X_bar(end-2);
Z_bar_3 = X_bar(end-1);
Z_bar_4 = X_bar(end);
X_bar_1 = X_bar(1);

%% Compute Derivatives
sigma_1 = k;
sigma_4 = k_0;
sigma_5 = alpha_1;
sigma_6 = alpha_2;

%% Compute PID Gains
K_P = -sigma_4;
K_I = sigma_1 * Z_bar_1 / (Z_bar_1 + Z_bar_2);
K_D = sigma_6 / sigma_5;
omega_c = eta * (Z_bar_1 + Z_bar_2);
omega_0 = sigma_5;
K_F = sigma_1 / (eta * (Z_bar_1 + Z_bar_2));
K_S = theta;

%% Construct Gains Structure
Gains.K_P = K_P;
Gains.K_I = K_I;
Gains.K_D = K_D;
Gains.omega_c = omega_c;
Gains.omega_0 = omega_0;
Gains.K_F = K_F;
Gains.K_S = K_S;

end

