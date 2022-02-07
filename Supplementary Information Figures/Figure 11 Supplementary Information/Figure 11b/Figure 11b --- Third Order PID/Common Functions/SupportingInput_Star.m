function [u, x] = SupportingInput_Star(Parameters, r)
%% Extract Parameters
gamma_F = Parameters.gamma_F;
kappa_F = Parameters.kappa_F;
gamma_1 = Parameters.gamma_1;
gamma_2 = Parameters.gamma_2;
gamma_3 = Parameters.gamma_3;
gamma_4 = Parameters.gamma_4;
gamma_5 = Parameters.gamma_5;
gamma_6 = Parameters.gamma_6;
k_1 = Parameters.k_1;
k_2 = Parameters.k_2;
k_3 = Parameters.k_3;
k_4 = Parameters.k_4;
k_5 = Parameters.k_5;

%% Compute Plant Fixed Point
x_6 = r;
x_5 = (gamma_6/k_5) * x_6;
x_4 = (gamma_5/k_4) * x_5;
x_3 = (gamma_4/k_3) * x_4;
x_2 = (gamma_3/k_2) * x_3;
x_1 = (x_2/k_1) * (gamma_2 + gamma_F*x_6 / (x_2 + kappa_F));
x = [x_1; x_2; x_3; x_4; x_5; x_6];

%% Compute Supporting Input
u = gamma_1 * x(1);

end

