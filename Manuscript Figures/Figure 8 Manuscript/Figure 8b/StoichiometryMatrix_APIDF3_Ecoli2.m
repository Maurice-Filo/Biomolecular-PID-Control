function S = StoichiometryMatrix_APIDF3_Ecoli2()
% Stoichiometry Matrix for Network APIDF3_Ecoli2
% 	 Species: 		 X = [X_1; X_2; Z_1; Z_2; Z_3]
% 	 Reactions: 	R1:         X_1					-->         phi				[gamma_1*X_1]
% 				    R2:         X_2					-->         phi				[gamma_2*X_2]
% 				    R3:         X_1					-->         X_1 +  X_2		[k_1*X_1]
% 				    R4:         phi					-->         Z_1				[mu]
% 				    R5:         phi					-->         Z_2				[theta*X_2]
% 				    R6:         Z_1 +  Z_2			-->         phi				[eta*Z_1*Z_2]
%                   R7:         phi                 -->         Z_3             [alpha_0 / (1 + (X_2/kappa_0)^n_0)]
%                   R8:         Z_3                 -->         phi             [gamma_0*Z_3]
% 				    R9:         phi					-->         X_1				[k*Z_1]
% 				    R10:        X_1					-->         phi				[(delta*X_2 + delta_0*Z_3) * (X_1/kappa) / (1 + (X_1/kappa)^n)]
% 				    R11:        Z_1					-->         phi				[delta_c*Z_1]
% 				    R12:        Z_2					-->         phi				[delta_c*Z_2]

S = [  -1,		0,		0,		0,		0,		0,		0,		0,		1,	   -1,		0,		0; ...
        0,	   -1,		1,		0,		0,		0,		0,		0,		0,		0,		0,		0; ...
		0,		0,		0,		1,		0,	   -1,		0,		0,		0,		0,	   -1,		0; ...
		0,		0,		0,		0,		1,	   -1,		0,		0,		0,		0,		0,	   -1; ...
		0,		0,		0,		0,		0,		0,		1,	   -1,		0,		0,		0,		0	];
end