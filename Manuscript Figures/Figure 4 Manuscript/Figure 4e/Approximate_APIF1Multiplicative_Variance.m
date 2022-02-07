function Variance = Approximate_APIF1Multiplicative_Variance(Parameters)
%% Extract Parameters
kappa = Parameters.kappa;
n = Parameters.n;
mu = Parameters.mu;
theta = Parameters.theta;
k = Parameters.k;
gamma_1 = Parameters.gamma_1;
gamma_2 = Parameters.gamma_2;
k_1 = Parameters.k_1;
r = mu/theta;

%% Compute Partial Derivatives
sigma_1 = k / (1 + (r/kappa)^n); 
sigma_3 = 0;
sigma_4 = (gamma_1*gamma_2/k_1) * n * (r/kappa)^n / (1 + (r/kappa)^n); 

%% Compute Variance
Variance = r * ( (gamma_1 + gamma_2 + sigma_3) * (gamma_1*gamma_2 + gamma_2*sigma_3 + sigma_1*k_1) + k_1*gamma_2 * (gamma_1 + sigma_4) ) / ...
            ( (gamma_1 + gamma_2 + sigma_3) * (gamma_1*gamma_2 + gamma_2*sigma_3 + k_1*sigma_4) - sigma_1 * k_1 * theta );
        
end

