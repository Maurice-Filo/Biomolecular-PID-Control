function X_bar = APIDF2_Deg_FixedPoint(Parameters, SupportingInput)
% u = k * z_1 - delta *  x_L * (x_1/kappa)^n / (1 + (x_1/kappa)^n)

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

%% Compute Supporting Input
[u, x] = SupportingInput(Parameters, r);

%% Compute Fixed Point
z_1 = (u + delta * r *  (x(1)/kappa)^n / (1 + (x(1)/kappa)^n ))/ k;
z_2 = (mu + beta*r) / (eta*z_1);

%% Stack Coordinates
X_bar = [x; z_1; z_2];
end

