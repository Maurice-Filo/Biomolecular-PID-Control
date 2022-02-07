function X_bar = APIDF4_Deg_FixedPoint(Parameters, SupportingInput)
% u = k * z_1 - (delta*x_L + alpha_1*z_3 + alpha_2*x_L) * (x_1/kappa)^n / (1 + (x_1/kappa)^n)

%% Extract Controller Parameters
mu = Parameters.mu;
theta = Parameters.theta;
k = Parameters.k;
delta = Parameters.delta;
mu_0 = Parameters.mu_0;
alpha_1 = Parameters.alpha_1;
alpha_2 = Parameters.alpha_2;
eta = Parameters.eta;
eta_0 = Parameters.eta_0;
kappa = Parameters.kappa;
n = Parameters.n;
delta_c = Parameters.delta_c;
r = mu / theta;

%% Compute Supporting Input
[u, x] = SupportingInput(Parameters, r);

%% Compute Fixed Point
z_3 = (mu_0 - alpha_2*r) / alpha_1;
z_4 = mu_0/(eta_0*z_3);
z_1 = (u + (delta*r + mu_0) *  (x(1)/kappa)^n / (1 + (x(1)/kappa)^n))/ k;
z_2 = mu / (eta*z_1);

%% Stack Coordinates
X_bar = [x; z_1; z_2; z_3; z_4];
end

