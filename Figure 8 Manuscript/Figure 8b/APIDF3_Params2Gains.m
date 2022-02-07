function Gains = APIDF3_Params2Gains(Parameters, SupportingInput)
% u = k * z_1   - (delta*x_L + delta_0*z_3) * (x_1/kappa)^n / (1 + (x_1/kappa)^n)

%% Extract Controller Parameters
mu = Parameters.mu;
theta = Parameters.theta;
k = Parameters.k;
eta = Parameters.eta;
delta = Parameters.delta;
delta_0 = Parameters.delta_0;
alpha_0 = Parameters.alpha_0;
gamma_0 = Parameters.gamma_0;
kappa_0 = Parameters.kappa_0;
kappa = Parameters.kappa;
n = Parameters.n;
n_0 = Parameters.n_0;
r = mu / theta;

%% Fixed Point
X_bar = APIDF3_IdealFixedPoint(Parameters, SupportingInput);
Z_bar_1 = X_bar(end-2);
Z_bar_2 = X_bar(end-1);
Z_bar_3 = X_bar(end);
X_bar_1 = X_bar(1);

%% Compute Derivatives
sigma_1 = k;
sigma_2 = delta_0 * (X_bar_1/kappa)^n / (1 + (X_bar_1/kappa)^n);
sigma_3 =  (delta*r + delta_0*Z_bar_3) * (n/X_bar_1) * (X_bar_1/kappa)^n / (1 + (X_bar_1/kappa)^n)^2;
sigma_4 = delta * (X_bar_1/kappa)^n / (1 + (X_bar_1/kappa)^n);
sigma_5 = 0;
sigma_6 = (n_0*alpha_0/r) * (r/kappa_0)^n_0 / (1 + (r/kappa_0)^n_0)^2;

%% Compute PID Gains
K_P = sigma_4 - sigma_2*sigma_6 / (gamma_0 + sigma_5);
K_I = sigma_1 * Z_bar_1 / (Z_bar_1 + Z_bar_2);
K_D = sigma_4 / (gamma_0 + sigma_5);
omega_c = eta * (Z_bar_1 + Z_bar_2);
omega_0 = gamma_0 + sigma_5;
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

