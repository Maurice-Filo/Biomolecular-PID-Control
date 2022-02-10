function S = StoichiometryMatrix_aPID2_Deg_FeedbackNet()
% Stoichiometry Matrix for Network aPID2_Deg_FeedbackNet
% 	 Species: 		 X = [X_1; X_2; X_3; X_4; X_5; X_6; Z_1; Z_2]
% 	 Reactions: 	R1:		 X_1					--> 	 X_1 +  X_2				[k_1*X_1]
% 				    R2:		 X_2					--> 	 X_2 +  X_3				[k_2*X_2]
% 				    R3:		 X_3					--> 	 X_3 +  X_4				[k_3*X_3]
% 				    R4:		 X_4					--> 	 X_4 +  X_5				[k_4*X_4]
% 				    R5:		 X_5					--> 	 X_5 +  X_6				[k_5*X_5]
% 				    R6:		 X_1					--> 	phi				[gamma_1*X_1]
% 				    R7:		 X_2					--> 	phi				[gamma_2*X_2 + gamma_F*X_6*X_2/(X_2 + kappa_F)]
% 				    R8:		 X_3					--> 	phi				[gamma_3*X_3]
% 				    R9:		 X_4					--> 	phi				[gamma_4*X_4]
% 				    R10:		 X_5					--> 	phi				[gamma_5*X_5]
% 				    R11:		 X_6					--> 	phi				[gamma_6*X_6]
% 				    R12:		phi					--> 	 Z_1				[mu + beta*X_6]
% 				    R13:		 X_6					--> 	 X_6 +  Z_2				[theta*X_6]
% 				    R14:		 Z_1 +  Z_2					--> 	phi				[eta*Z_1*Z_2]
% 				    R15:		phi					--> 	 X_1				[k*Z_1]
% 				    R16:		 X_1					--> 	phi				[delta*X_6*X_1/(X_1 + kappa_1)]
% 				    R17:		 Z_1					--> 	phi				[gamma_c*Z_1]
% 				    R18:		 Z_2					--> 	phi				[gamma_c*Z_2]
S = [	0,		0,		0,		0,		0,		-1,		0,		0,		0,		0,		0,		0,		0,		0,		1,		-1,		0,		0; ...
		1,		0,		0,		0,		0,		0,		-1,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0; ...
		0,		1,		0,		0,		0,		0,		0,		-1,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0; ...
		0,		0,		1,		0,		0,		0,		0,		0,		-1,		0,		0,		0,		0,		0,		0,		0,		0,		0; ...
		0,		0,		0,		1,		0,		0,		0,		0,		0,		-1,		0,		0,		0,		0,		0,		0,		0,		0; ...
		0,		0,		0,		0,		1,		0,		0,		0,		0,		0,		-1,		0,		0,		0,		0,		0,		0,		0; ...
		0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		1,		0,		-1,		0,		0,		-1,		0; ...
		0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		0,		1,		-1,		0,		0,		0,		-1	];
end