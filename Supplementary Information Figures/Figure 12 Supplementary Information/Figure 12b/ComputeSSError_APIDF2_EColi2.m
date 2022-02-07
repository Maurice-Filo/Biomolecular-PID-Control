function [e, Output] = ComputeSSError_APIDF2_EColi2(Parameters, r)
%% Extract Parameters
gamma_1 = Parameters.gamma_1;
gamma_2 = Parameters.gamma_2;
k_1 = Parameters.k_1;
mu = Parameters.mu;
theta = Parameters.theta;
eta = Parameters.eta;
k = Parameters.k;
delta_c = Parameters.delta_c;
delta = Parameters.delta;
kappa = Parameters.kappa;
beta = Parameters.beta;

%% Computing the Setpoint
Output = (k*k_1*((gamma_1*gamma_2*delta_c^4 + 2*gamma_1*gamma_2*delta_c^2*eta*mu + 4*k*k_1*theta*delta_c*eta*mu + gamma_1*gamma_2*eta^2*mu^2)^(1/2) - delta_c^2*gamma_1^(1/2)*gamma_2^(1/2) + eta*gamma_1^(1/2)*gamma_2^(1/2)*mu))/(eta*gamma_1^(1/2)*gamma_2^(1/2)*(2*delta_c*gamma_1*gamma_2 + 2*k*k_1*theta));
e = Output - r;
end

