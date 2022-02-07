function a = PropensityFunction_APIDF2_Deg_GeneExp(x, Parameters)
% Propensity for Network APIDF2_Deg_GeneExp
% 	 Species: 		 X = [X_1; X_2; Z_1; Z_2]
% 	 Reactions: 	R1:		 X_1					--> 	phi				[gamma_1*X_1]
% 				    R2:		 X_2					--> 	phi				[gamma_2*X_2]
% 				    R3:		 X_1					--> 	 X_1 +  X_2				[k_1*X_1]
% 				    R4:		phi					--> 	 Z_1				[mu + beta*X_2]
% 				    R5:		phi					--> 	 Z_2				[theta*X_2]
% 				    R6:		 Z_1 +  Z_2					--> 	phi				[eta*Z_1*Z_2]
% 				    R7:		phi					--> 	 X_1				[k*Z_1]
% 				    R8:		 X_1					--> 	phi				[delta*X_2 * (X_1/kappa)^n / (1 + (X_1/kappa)^n)]
% 				    R9:		 Z_1					--> 	phi				[delta_c*Z_1]
% 				    R10:		 Z_2					--> 	phi				[delta_c*Z_2]

%% Extract Parameters
gamma_1 = Parameters.gamma_1;
gamma_2 = Parameters.gamma_2;
k_1 = Parameters.k_1;
mu = Parameters.mu;
theta = Parameters.theta;
beta = Parameters.beta;
k = Parameters.k;
eta = Parameters.eta;
delta = Parameters.delta;
kappa = Parameters.kappa;
n = Parameters.n;
delta_c = Parameters.delta_c;

%% Extract State Variables
X_1 = x(1);
X_2 = x(2);
Z_1 = x(3);
Z_2 = x(4);

%% Propensities
a = [ ...
		gamma_1*X_1; ...
		gamma_2*X_2; ...
		k_1*X_1; ...
		mu + beta*X_2; ...
		theta*X_2; ...
		eta*Z_1*Z_2; ...
		k*Z_1; ...
		delta*X_2 * (X_1/kappa)^n / (1 + (X_1/kappa)^n); ...
		delta_c*Z_1; ...
		delta_c*Z_2 ...
	];
end