function S = StoichiometryMatrix_OL_FeedbackNetwork()
% Stoichiometry Matrix for Network OL_FeedbackNetwork
% 	 Species: 		 X = [X_1; X_2; X_3; X_4; X_5; X_6]
% 	 Reactions: 	R1:		phi					--> 	 X_1				[k_0]
% 				    R2:		 X_1					--> 	 X_1 +  X_2				[k_1*X_1]
% 				    R3:		 X_2					--> 	 X_2 +  X_3				[k_2*X_2]
% 				    R4:		 X_3					--> 	 X_3 +  X_4				[k_3*X_3]
% 				    R5:		 X_4					--> 	 X_4 +  X_5				[k_4*X_4]
% 				    R6:		 X_5					--> 	 X_5 +  X_6				[k_5*X_5]
% 				    R7:		 X_1					--> 	phi				[gamma_1*X_1]
% 				    R8:		 X_2					--> 	phi				[gamma_2*X_2 + gamma_F*X_6*X_2/(X_2 + kappa_F)]
% 				    R9:		 X_3					--> 	phi				[gamma_3*X_3]
% 				    R10:		 X_4					--> 	phi				[gamma_4*X_4]
% 				    R11:		 X_5					--> 	phi				[gamma_5*X_5]
% 				    R12:		 X_6					--> 	phi				[gamma_6*X_6]
S = [	1,		0,		0,		0,		0,		0,		-1,		0,		0,		0,		0,		0; ...
		0,		1,		0,		0,		0,		0,		0,		-1,		0,		0,		0,		0; ...
		0,		0,		1,		0,		0,		0,		0,		0,		-1,		0,		0,		0; ...
		0,		0,		0,		1,		0,		0,		0,		0,		0,		-1,		0,		0; ...
		0,		0,		0,		0,		1,		0,		0,		0,		0,		0,		-1,		0; ...
		0,		0,		0,		0,		0,		1,		0,		0,		0,		0,		0,		-1	];
end