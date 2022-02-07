function [A, B, W, C] = APIDF4H_Deg_Jacobians(Parameters, FixedPoint, SupportingInput, JacobiansPlant)
% u = k * z_1 + k_0 * x_L - (alpha_1*Z_3 + alpha_2*X_6) * (x_1/kappa)^n / (1 + (x_1/kappa)^n)

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

%% Compute Open Loop Jacobian
[A_OL, ~, W_OL, C_OL] = JacobiansPlant(Parameters, r);

%% Compute Closed Loop Fixed Point
X_bar = FixedPoint(Parameters, SupportingInput);
X_bar_1 = X_bar(1);
Z_bar_1 = X_bar(end-3);
Z_bar_2 = X_bar(end-2);
Z_bar_3 = X_bar(end-1);
Z_bar_4 = X_bar(end);

%% Compute Actuation Derivatives
T = (X_bar_1/kappa)^n / (1 + (X_bar_1/kappa)^n);
sigma_1 = k;
sigma_3 = mu_0 * (n/X_bar_1) * (X_bar_1/kappa)^n / (1 + (X_bar_1/kappa)^n)^2;
sigma_4 = k_0;
sigma_5 = alpha_1;
sigma_6 = alpha_2;

%% Compute Closed Loop Jacobians
L = size(A_OL,1);
A = [A_OL,             	    zeros(L,4); ...
   	 zeros(1,L),           -eta*Z_bar_2,           -eta*Z_bar_1,            0,                      0; ...
     zeros(1,L),           -eta*Z_bar_2,           -eta*Z_bar_1,            0,                      0; ...
     zeros(1,L),            0,                      0,                      -eta_0*Z_bar_4,         -eta_0*Z_bar_3; ...
     zeros(1,L),            0,                      0,                      -eta_0*Z_bar_4,         -eta_0*Z_bar_3];
A(1,1) = A(1,1) - sigma_3;
A(1,L) = A(1,L) + sigma_4 - sigma_6*T;
A(1,L+1) = sigma_1;
A(1,L+3) = -sigma_5*T;
A(L+2,L) = theta;
A(L+4,L) = sigma_6;
A(L+4,L+3) = A(L+4,L+3) + sigma_5;
B = [0; 0; 1; 0; 0; 0];
W = [W_OL; 0; 0; 0; 0];
C = [C_OL, 0, 0, 0, 0];
end

