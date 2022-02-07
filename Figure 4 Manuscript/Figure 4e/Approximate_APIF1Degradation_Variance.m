function Variance = Approximate_APIF1Degradation_Variance(Parameters)
%% Extract Parameters
delta = Parameters.delta;
kappa_1 = Parameters.kappa_1;
n = Parameters.n;
mu = Parameters.mu;
theta = Parameters.theta;
k = Parameters.k;
gamma_1 = Parameters.gamma_1;
gamma_2 = Parameters.gamma_2;
k_1 = Parameters.k_1;
r = mu/theta;

%% Compute Partial Derivatives
x_bar_1 = gamma_2*r/k_1;
sigma_1 = k;
sigma_3 = delta*r^n * kappa_1 / (x_bar_1 + kappa_1)^2;
sigma_4 = n * delta * r^(n-1)  * x_bar_1 / (x_bar_1 + kappa_1); 
h_minus = delta*r^n * (x_bar_1) / (kappa_1 + x_bar_1);

%% Compute Variance
Variance = r * ( (gamma_1 + gamma_2 + sigma_3) * (gamma_1*gamma_2 + gamma_2*sigma_3 + sigma_1*k_1) + k_1*gamma_2 * (gamma_1 + sigma_4) + (k_1^2/r) * h_minus ) / ...
            ( (gamma_1 + gamma_2 + sigma_3) * (gamma_1*gamma_2 + gamma_2*sigma_3 + k_1*sigma_4) - sigma_1 * k_1 * theta );
end

