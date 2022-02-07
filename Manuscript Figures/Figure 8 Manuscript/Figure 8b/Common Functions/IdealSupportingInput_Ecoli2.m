function [u, x] = IdealSupportingInput_Ecoli2(Parameters, r)
%% Extract Plant Parameters
gamma_1 = Parameters.gamma_1;
gamma_2 = Parameters.gamma_2;
k_1 = Parameters.k_1;

%% Compute Plant Fixed Point
x = zeros(2,1);
x(2) = r;
x(1) = gamma_2 * x(2) / k_1;

%% Compute Supporting Input
u = gamma_1 * x(1);

end

