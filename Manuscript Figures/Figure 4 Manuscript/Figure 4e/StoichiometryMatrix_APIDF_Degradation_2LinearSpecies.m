function S = StoichiometryMatrix_APIDF_Degradation_2LinearSpecies()
% Stoichiometry Matrix for Network APIDF_Degradation_2LinearSpecies
% 	 Species: 		 X = [X_1; X_2; Z_1; Z_2]
% 	 Reactions: 	R1:		 X_1					--> 	 X_1 +  X_2				[k_1*X_1]
% 				    R2:		 X_1					--> 	phi				[gamma_1*X_1]
% 				    R3:		 X_2					--> 	phi				[gamma_2*X_2]
% 				    R4:		phi					--> 	 Z_1				[mu + beta*X_2]
% 				    R5:		 X_2					--> 	 X_2 +  Z_2				[theta*X_2]
% 				    R6:		 Z_1 +  Z_2					--> 	phi				[eta*Z_1*Z_2]
% 				    R7:		 Z_1					--> 	 X_1 +  Z_1				[k*Z_1]
% 				    R8:		 X_1 +  X_2					--> 	 X_2				[delta*X_1*X_2^n]
% 				    R9:		 Z_1					--> 	phi				[gamma_c*Z_1]
% 				    R10:		 Z_2					--> 	phi				[gamma_c*Z_2]
S = [	0,		-1,		0,		0,		0,		0,		1,		-1,		0,		0; ...
		1,		0,		-1,		0,		0,		0,		0,		0,		0,		0; ...
		0,		0,		0,		1,		0,		-1,		0,		0,		-1,		0; ...
		0,		0,		0,		0,		1,		-1,		0,		0,		0,		-1	];
end