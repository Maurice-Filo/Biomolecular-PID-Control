function a = PropensityFunction_OL_Star(x, Parameters)
% Propensity for Network OL_Star
% 	 Species: 		 X = [X_1; X_2; X_3; X_4; X_5; X_6]
% 	 Reactions: 	R1:		phi					--> 	 X_1				[k]
% 				    R2:		 X_1					--> 	phi				[gamma_1*X_1]
% 				    R3:		 X_2					--> 	phi				[gamma_2*X_2]
% 				    R4:		 X_3					--> 	phi				[gamma_3*X_3]
% 				    R5:		 X_4					--> 	phi				[gamma_4*X_4]
% 				    R6:		 X_5					--> 	phi				[gamma_5*X_5]
% 				    R7:		 X_6					--> 	phi				[gamma_6*X_6]
% 				    R8:		 X_1					--> 	 X_1 +  X_2				[k_1*X_1]
% 				    R9:		 X_2					--> 	 X_2 +  X_3				[k_2*X_2]
% 				    R10:		 X_3					--> 	 X_3 +  X_4				[k_3*X_3]
% 				    R11:		 X_4					--> 	 X_4 +  X_5				[k_4*X_4]
% 				    R12:		 X_5					--> 	 X_5 +  X_6				[k_5*X_5]
% 				    R13:		 X_2					--> 	phi				[gamma_F*X_2*X_6/(X_2 + kappa_F)]

%% Extract Parameters
k = Parameters.k;
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

%% Extract State Variables
X_1 = x(1);
X_2 = x(2);
X_3 = x(3);
X_4 = x(4);
X_5 = x(5);
X_6 = x(6);

%% Propensities
a = [ ...
		k; ...
		gamma_1*X_1; ...
		gamma_2*X_2; ...
		gamma_3*X_3; ...
		gamma_4*X_4; ...
		gamma_5*X_5; ...
		gamma_6*X_6; ...
		k_1*X_1; ...
		k_2*X_2; ...
		k_3*X_3; ...
		k_4*X_4; ...
		k_5*X_5; ...
		gamma_F*X_2*X_6/(X_2 + kappa_F) ...
	];
end