function a = PropensityFunction_APICType1Additive_2SpeciesLinear(x, Parameters)
% Propensity for Network APICType1Additive_2SpeciesLinear
% 	 Species: 		 X = [X_1; X_2; Z_1; Z_2]
% 	 Reactions: 	R1:		 X_1					--> 	 X_1 +  X_2				[k_1*X_1]
% 				    R2:		 X_1					--> 	phi				[gamma_1*X_1]
% 				    R3:		 X_2					--> 	phi				[gamma_2*X_2]
% 				    R4:		phi					--> 	 Z_1				[mu]
% 				    R5:		 X_2					--> 	 X_2 +  Z_2				[theta*X_2]
% 				    R6:		 Z_1 +  Z_2					--> 	phi				[eta*Z_1*Z_2]
% 				    R7:		phi					--> 	 X_1				[k*Z_1 + alpha/(1+(X_2/kappa)^n)]
% 				    R8:		 Z_1					--> 	phi				[delta*Z_1]
% 				    R9:		 Z_2					--> 	phi				[delta*Z_2]

%% Extract Parameters
k_1 = Parameters.k_1;
gamma_1 = Parameters.gamma_1;
gamma_2 = Parameters.gamma_2;
mu = Parameters.mu;
theta = Parameters.theta;
eta = Parameters.eta;
k = Parameters.k;
alpha = Parameters.alpha;
kappa = Parameters.kappa;
n = Parameters.n;
delta = Parameters.delta;

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
		mu; ...
		theta*X_2; ...
		eta*Z_1*Z_2; ...
		k*Z_1 + alpha/(1+(X_2/kappa)^n); ...
		delta*Z_1; ...
		delta*Z_2 ...
	];
end