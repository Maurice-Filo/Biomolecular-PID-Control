function a = PropensityFunction_aPID3_Deg_GeneExp(x, Parameters)
% Propensity for Network aPID3_Deg_GeneExp
% 	 Species: 		 X = [X_1; X_2; Z_1; Z_2; Z_3]
% 	 Reactions: 	R1:		 X_1					--> 	 X_1 +  X_2				[k_1*X_1]
% 				    R2:		 X_1					--> 	phi				[gamma_1*X_1]
% 				    R3:		 X_2					--> 	phi				[gamma_2*X_2]
% 				    R4:		phi					--> 	 Z_1				[mu]
% 				    R5:		 X_2					--> 	 X_2 +  Z_2				[theta*X_2]
% 				    R6:		 Z_1 +  Z_2					--> 	phi				[eta*Z_1*Z_2]
% 				    R7:		phi					--> 	 Z_3				[alpha_0 / (1 + (X_2/kappa_0)^n_0)]
% 				    R8:		 Z_3					--> 	phi				[gamma_0*Z_3]
% 				    R9:		phi					--> 	 X_1				[k*Z_1]
% 				    R10:		 X_1					--> 	phi				[(delta*X_2 + delta_0*Z_3) * X_1 / (X_1 + kappa_1)]
% 				    R11:		 Z_1					--> 	phi				[gamma_c*Z_1]
% 				    R12:		 Z_2					--> 	phi				[gamma_c*Z_2]

%% Extract Parameters
k_1 = Parameters.k_1;
gamma_1 = Parameters.gamma_1;
gamma_2 = Parameters.gamma_2;
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
Z_1 = x(3);
Z_2 = x(4);
Z_3 = x(5);

%% Propensities
a = [ ...
		k_1*X_1; ...
		gamma_1*X_1; ...
		gamma_2*X_2; ...
		mu; ...
		theta*X_2; ...
		eta*Z_1*Z_2; ...
		alpha_0 / (1 + (X_2/kappa_0)^n_0); ...
		gamma_0*Z_3; ...
		k*Z_1; ...
		(delta*X_2 + delta_0*Z_3) * X_1 / (X_1 + kappa_1); ...
		gamma_c*Z_1; ...
		gamma_c*Z_2 ...
	];
end