function a = PropensityFunction_OL(x, Parameters)
% Propensity Function for Open Loop System
%   Species:        X = [X_1; X_2]
%   Reactions:      R1:     phi     --> X_1         [k]
%                   R2:     X_1     --> phi         [gamma_1]
%                   R3:     X_1     --> X_1 + X_2   [k_1]
%                   R4:     X_2     --> phi         [gamma_2]

    a = [Parameters.k; ...
         Parameters.gamma_1 * x(1); ...
         Parameters.k_1 * x(1); ...
         Parameters.gamma_2 * x(2)];
end