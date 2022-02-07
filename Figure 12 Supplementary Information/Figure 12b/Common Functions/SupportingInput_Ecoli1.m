function [u, x] = SupportingInput_Ecoli1(Parameters, r)
%% Extract Plant Parameters
gamma = Parameters.gamma;

%% Compute Plant Fixed Point
x = r;

%% Compute Supporting Input
u = gamma * x;

end

