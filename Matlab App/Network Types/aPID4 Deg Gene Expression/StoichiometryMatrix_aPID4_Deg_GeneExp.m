function S = StoichiometryMatrix_aPID4_Deg_GeneExp()
% Stoichiometry Matrix for Network aPID4_Deg_GeneExp
% 	 Species: 		 X = [X_1; X_2; Z_1; Z_2; Z_3; Z_4]
% 	 Reactions: 	R1:		 X_1					--> 	 X_1 +  X_2				[k_1*X_1]
% 				    R2:		 X_1					--> 	phi				[gamma_1*X_1]
% 				    R3:		 X_2					--> 	phi				[gamma_2*X_2]
% 				    R4:		phi					--> 	 Z_1				[mu]
% 				    R5:		 X_2					--> 	 X_2 +  Z_2				[theta*X_2]
% 				    R6:		 Z_1 +  Z_2					--> 	phi				[eta*Z_1*Z_2]
% 				    R7:		phi					--> 	 Z_3				[mu_0]
% 				    R8:		phi					--> 	 Z_4				[theta_0 * (alpha_1*Z_3 + alpha_2*X_2)]
% 				    R9:		 Z_3 +  Z_4					--> 	phi				[eta_0*Z_3*Z_4]
% 				    R10:		phi					--> 	 X_1				[k*Z_1]
% 				    R11:		 X_1					--> 	phi				[(delta*X_2 + alpha_1*Z_3 + alpha_2*X_2) * X_1 / (X_1 + kappa_1)]
% 				    R12:		 Z_1					--> 	phi				[gamma_c*Z_1]
% 				    R13:		 Z_2					--> 	phi				[gamma_c*Z_2]
% 				    R14:		 Z_3					--> 	phi				[gamma_c*Z_3]
% 				    R15:		 Z_4					--> 	phi				[gamma_c*Z_4]
S = [	0,		-1,		0,		0,		0,		0,		0,		0,		0,		1,		-1,		0,		0,		0,		0; ...
		1,		0,		-1,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0; ...
		0,		0,		0,		1,		0,		-1,		0,		0,		0,		0,		0,		-1,		0,		0,		0; ...
		0,		0,		0,		0,		1,		-1,		0,		0,		0,		0,		0,		0,		-1,		0,		0; ...
		0,		0,		0,		0,		0,		0,		1,		0,		-1,		0,		0,		0,		0,		-1,		0; ...
		0,		0,		0,		0,		0,		0,		0,		1,		-1,		0,		0,		0,		0,		0,		-1	];
end