function [A, B, W, C] = Jacobians_Star(Parameters, r)
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

%% Compute Fixed Point
[~, x] = SupportingInput_Star(Parameters, r);
X_bar_1 = x(1);
X_bar_2 = x(2);
X_bar_6 = r;

%% Compute Plant Jacobians
A = [-gamma_1,                                                                                                 0,        0,        0,        0,                                      0; ...
          k_1, (X_bar_2*X_bar_6*gamma_F)/(X_bar_2 + kappa_F)^2 - (X_bar_6*gamma_F)/(X_bar_2 + kappa_F) - gamma_2,        0,        0,        0, -(X_bar_2*gamma_F)/(X_bar_2 + kappa_F); ...
            0,                                                                                               k_2, -gamma_3,        0,        0,                                      0; ...
            0,                                                                                                 0,      k_3, -gamma_4,        0,                                      0; ...
            0,                                                                                                 0,        0,      k_4, -gamma_5,                                      0; ...
            0,                                                                                                 0,        0,        0,      k_5,                               -gamma_6];
B = [1; 0; 0; 0; 0; 0];
W = [0; X_bar_1; 0; 0; 0; 0];
C = [0, 0, 0, 0, 0, 1];
end

