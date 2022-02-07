function X_bar = APIDF3_Deg_FixedPoint(Parameters, SupportingInput)
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

%% Compute Supporting Input
[u, x] = SupportingInput(Parameters, r);

%% Compute Fixed Point
z_3 = (alpha_0/gamma_0) / (1 + (r/kappa_0)^n_0);
z_1 = (u + (delta*r + delta_0*z_3) *  (x(1)/kappa)^n / (1 + (x(1)/kappa)^n))/ k;
z_2 = mu / (eta*z_1);

%% Stack Coordinates
X_bar = [x; z_1; z_2; z_3];
end

