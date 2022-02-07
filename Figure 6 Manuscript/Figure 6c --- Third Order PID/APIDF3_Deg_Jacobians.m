function [A, B, W, C] = APIDF3_Deg_Jacobians(Parameters, FixedPoint, SupportingInput, JacobiansPlant)
% u = k * z_1 - (delta*x_L + delta_0*z_3) * (x_1/kappa)^n / (1 + (x_1/kappa)^n)

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

%% Compute Open Loop Jacobian
[A_OL, ~, W_OL, C_OL] = JacobiansPlant(Parameters, r);

%% Compute Closed Loop Fixed Point
X_bar = FixedPoint(Parameters, SupportingInput);
X_bar_1 = X_bar(1);
Z_bar_1 = X_bar(end-2);
Z_bar_2 = X_bar(end-1);
Z_bar_3 = X_bar(end);

%% Compute Actuation Derivatives
T = (X_bar_1/kappa)^n / (1 + (X_bar_1/kappa)^n);
sigma_1 = k;
sigma_2 = delta_0 * T;
sigma_3 = (delta*r + delta_0*Z_bar_3) * (n/X_bar_1) * (X_bar_1/kappa)^n / (1 + (X_bar_1/kappa)^n)^2;
sigma_4 = delta * T;
sigma_5 = 0;
sigma_6 = (n_0*alpha_0/r) * (r/kappa_0)^n_0 / (1 + (r/kappa_0)^n_0)^2;

%% Compute Closed Loop Jacobians
L = size(A_OL,1);
A = [A_OL,             	    zeros(L,3); ...
   	 zeros(1,L),           -eta*Z_bar_2,           -eta*Z_bar_1,            0; ...
     zeros(1,L),           -eta*Z_bar_2,           -eta*Z_bar_1,            0; ...
     zeros(1,L),            0,                      0,                      -gamma_0];
A(1,1) = A(1,1) - sigma_3;
A(1,L) = A(1,L) - sigma_4;
A(1,L+1) = sigma_1;
A(1,L+3) = -sigma_2;
A(L+2,L) = theta;
A(L+3,L) = -sigma_6;
A(L+3,L+3) = -sigma_5 - gamma_0;
B = [0; 0; 1; 0; 0];
W = [W_OL; 0; 0; 0];
C = [C_OL, 0, 0, 0];

end

