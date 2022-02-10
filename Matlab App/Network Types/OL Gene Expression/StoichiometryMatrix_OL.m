function S = StoichiometryMatrix_OL()
% Stoichiometry Matrix for Open Loop System
%   Species:        X = [M; P]
%   Reactions:      R1:     phi --> M       [k_M]
%                   R2:     M   --> phi     [gamma_M]
%                   R3:     M   --> M + P   [k_P]
%                   R4:     P   --> phi     [gamma_P]
S = [1,       -1,     0,      0; ...
     0,       0,      1,      -1];
end

