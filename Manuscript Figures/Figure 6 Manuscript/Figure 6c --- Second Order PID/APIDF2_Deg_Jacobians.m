function [A, B, W, C] = APIDF2_Deg_Jacobians(Parameters, FixedPoint, SupportingInput, JacobiansPlant)
% u = k * z_1 - delta*x_L * (x_1/kappa)^n / (1 + (x_1/kappa)^n)

%% Extract Controller Parameters
mu = Parameters.mu;
theta = Parameters.theta;
beta = Parameters.beta;
k = Parameters.k;
eta = Parameters.eta;
kappa = Parameters.kappa;
delta = Parameters.delta;
n = Parameters.n;
r = mu / (theta - beta);

%% Compute Open Loop Jacobian
[A_OL, ~, W_OL, C_OL] = JacobiansPlant(Parameters, r);

%% Compute Closed Loop Fixed Point
X_bar = FixedPoint(Parameters, SupportingInput);
X_bar_1 = X_bar(1);
Z_bar_1 = X_bar(end-1);
Z_bar_2 = X_bar(end);

%% Compute Actuation Derivatives
T = (X_bar_1/kappa)^n / (1 + (X_bar_1/kappa)^n);
sigma_1 = k;
sigma_3 = (delta*r) * (n/X_bar_1) * (X_bar_1/kappa)^n / (1 + (X_bar_1/kappa)^n)^2;
sigma_4 = delta * T;

%% Compute Closed Loop Jacobians
L = size(A_OL,1);
A = [A_OL,             	    zeros(L,2); ...
   	 zeros(1,L),           -eta*Z_bar_2,           -eta*Z_bar_1; ...
     zeros(1,L),           -eta*Z_bar_2,           -eta*Z_bar_1];
A(1,1) = A(1,1) - sigma_3;
A(1,L) = A(1,L) - sigma_4;
A(L+1,L) = beta;
A(1,L+1) = sigma_1;
A(L+2,L) = theta;
B = [0; 0; 1; 0];
W = [W_OL; 0; 0];
C = [C_OL, 0, 0];

end

