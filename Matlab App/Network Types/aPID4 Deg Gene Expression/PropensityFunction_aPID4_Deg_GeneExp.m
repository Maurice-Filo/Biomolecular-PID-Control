function a = PropensityFunction_aPID4_Deg_GeneExp(x, Parameters)
% Propensity for Network aPID4_Deg_GeneExp
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

%% Extract Parameters
k_1 = Parameters.k_1;
gamma_1 = Parameters.gamma_1;
gamma_2 = Parameters.gamma_2;
mu = Parameters.mu;
theta = Parameters.theta;
eta = Parameters.eta;
k = Parameters.k;
delta = Parameters.delta;
alpha_1 = Parameters.alpha_1;
alpha_2 = Parameters.alpha_2;
kappa_1 = Parameters.kappa_1;
mu_0 = Parameters.mu_0;
theta_0 = Parameters.theta_0;
eta_0 = Parameters.eta_0;
gamma_c = Parameters.gamma_c;

%% Extract State Variables
X_1 = x(1);
X_2 = x(2);
Z_1 = x(3);
Z_2 = x(4);
Z_3 = x(5);
Z_4 = x(6);

%% Propensities
a = [ ...
		k_1*X_1; ...
		gamma_1*X_1; ...
		gamma_2*X_2; ...
		mu; ...
		theta*X_2; ...
		eta*Z_1*Z_2; ...
		mu_0; ...
		theta_0 * (alpha_1*Z_3 + alpha_2*X_2); ...
		eta_0*Z_3*Z_4; ...
		k*Z_1; ...
		(delta*X_2 + alpha_1*Z_3 + alpha_2*X_2) * X_1 / (X_1 + kappa_1); ...
		gamma_c*Z_1; ...
		gamma_c*Z_2; ...
		gamma_c*Z_3; ...
		gamma_c*Z_4 ...
	];
end