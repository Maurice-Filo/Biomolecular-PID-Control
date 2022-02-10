function a = PropensityFunction_aPID2_Deg_GeneExp(x, Parameters)
% Propensity for Network aPID2_Deg_GeneExp
% 	 Species: 		 X = [X_1; X_2; Z_1; Z_2]
% 	 Reactions: 	R1:		 X_1					--> 	 X_1 +  X_2				[k_1*X_1]
% 				    R2:		 X_1					--> 	phi				[gamma_1*X_1]
% 				    R3:		 X_2					--> 	phi				[gamma_2*X_2]
% 				    R4:		phi					--> 	 Z_1				[mu + beta*X_2]
% 				    R5:		 X_2					--> 	 X_2 +  Z_2				[theta*X_2]
% 				    R6:		 Z_1 +  Z_2					--> 	phi				[eta*Z_1*Z_2]
% 				    R7:		phi					--> 	 X_1				[k*Z_1]
% 				    R8:		 X_1					--> 	phi				[delta*X_2*X_1/(X_1 + kappa_1)]
% 				    R9:		 Z_1					--> 	phi				[gamma_c*Z_1]
% 				    R10:		 Z_2					--> 	phi				[gamma_c*Z_2]

%% Extract Parameters
k_1 = Parameters.k_1;
gamma_1 = Parameters.gamma_1;
gamma_2 = Parameters.gamma_2;
mu = Parameters.mu;
theta = Parameters.theta;
beta = Parameters.beta;
eta = Parameters.eta;
k = Parameters.k;
delta = Parameters.delta;
kappa_1 = Parameters.kappa_1;
gamma_c = Parameters.gamma_c;

%% Extract State Variables
X_1 = x(1);
X_2 = x(2);
Z_1 = x(3);
Z_2 = x(4);

%% Propensities
a = [ ...
		k_1*X_1; ...
		gamma_1*X_1; ...
		gamma_2*X_2; ...
		mu + beta*X_2; ...
		theta*X_2; ...
		eta*Z_1*Z_2; ...
		k*Z_1; ...
		delta*X_2*X_1/(X_1 + kappa_1); ...
		gamma_c*Z_1; ...
		gamma_c*Z_2 ...
	];
end