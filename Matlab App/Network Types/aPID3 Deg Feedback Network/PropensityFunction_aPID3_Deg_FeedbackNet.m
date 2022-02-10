function a = PropensityFunction_aPID3_Deg_FeedbackNet(x, Parameters)
% Propensity for Network aPID3_Deg_FeedbackNet
% 	 Species: 		 X = [X_1; X_2; X_3; X_4; X_5; X_6; Z_1; Z_2; Z_3]
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
% 				    R12:		phi					--> 	 Z_1				[mu]
% 				    R13:		 X_6					--> 	 X_6 +  Z_2				[theta*X_6]
% 				    R14:		 Z_1 +  Z_2					--> 	phi				[eta*Z_1*Z_2]
% 				    R15:		phi					--> 	 Z_3				[alpha_0 / (1 + (X_6/kappa_0)^n_0)]
% 				    R16:		 Z_3					--> 	phi				[gamma_0*Z_3]
% 				    R17:		phi					--> 	 X_1				[k*Z_1]
% 				    R18:		 X_1					--> 	phi				[(delta*X_6 + delta_0*Z_3)*X_1/(X_1 + kappa_1)]
% 				    R19:		 Z_1					--> 	phi				[gamma_c*Z_1]
% 				    R20:		 Z_2					--> 	phi				[gamma_c*Z_2]

%% Extract Parameters
k_1 = Parameters.k_1;
k_2 = Parameters.k_2;
k_3 = Parameters.k_3;
k_4 = Parameters.k_4;
k_5 = Parameters.k_5;
gamma_1 = Parameters.gamma_1;
gamma_2 = Parameters.gamma_2;
gamma_3 = Parameters.gamma_3;
gamma_4 = Parameters.gamma_4;
gamma_5 = Parameters.gamma_5;
gamma_6 = Parameters.gamma_6;
gamma_F = Parameters.gamma_F;
kappa_F = Parameters.kappa_F;
mu = Parameters.mu;
theta = Parameters.theta;
eta = Parameters.eta;
k = Parameters.k;
delta = Parameters.delta;
delta_0 = Parameters.delta_0;
kappa_1 = Parameters.kappa_1;
alpha_0 = Parameters.alpha_0;
kappa_0 = Parameters.kappa_0;
n_0 = Parameters.n_0;
gamma_0 = Parameters.gamma_0;
gamma_c = Parameters.gamma_c;

%% Extract State Variables
X_1 = x(1);
X_2 = x(2);
X_3 = x(3);
X_4 = x(4);
X_5 = x(5);
X_6 = x(6);
Z_1 = x(7);
Z_2 = x(8);
Z_3 = x(9);

%% Propensities
a = [ ...
		k_1*X_1; ...
		k_2*X_2; ...
		k_3*X_3; ...
		k_4*X_4; ...
		k_5*X_5; ...
		gamma_1*X_1; ...
		gamma_2*X_2 + gamma_F*X_6*X_2/(X_2 + kappa_F); ...
		gamma_3*X_3; ...
		gamma_4*X_4; ...
		gamma_5*X_5; ...
		gamma_6*X_6; ...
		mu; ...
		theta*X_6; ...
		eta*Z_1*Z_2; ...
		alpha_0 / (1 + (X_6/kappa_0)^n_0); ...
		gamma_0*Z_3; ...
		k*Z_1; ...
		(delta*X_6 + delta_0*Z_3)*X_1/(X_1 + kappa_1); ...
		gamma_c*Z_1; ...
		gamma_c*Z_2 ...
	];
end