function a = PropensityFunction_APIDF4_Deg_GeneExp(x, Parameters)
% Propensity for Network APIDF4_Deg_GeneExp
% 	 Species: 		 X = [X_1; X_2; Z_1; Z_2; Z_3; Z_4]
% 	 Reactions: 	R1:		 X_1					--> 	phi				[gamma_1*X_1]
% 				    R2:		 X_2					--> 	phi				[gamma_2*X_2]
% 				    R3:		 X_1					--> 	 X_1 +  X_2				[k_1*X_1]
% 				    R4:		phi					--> 	 Z_1				[mu]
% 				    R5:		phi					--> 	 Z_2				[theta*X_2]
% 				    R6:		 Z_1 +  Z_2					--> 	phi				[eta*Z_1*Z_2]
% 				    R7:		phi					--> 	 Z_3				[mu_0]
% 				    R8:		phi					--> 	 Z_4				[alpha_1*Z_3 + alpha_2*X_2]
% 				    R9:		 Z_3 +  Z_4					--> 	phi				[eta_0*Z_3*Z_4]
% 				    R10:		phi					--> 	 X_1				[k*Z_1]
% 				    R11:		 X_1					--> 	phi				[(delta*X_2 + alpha_1*Z_3 + alpha_2*X_2) * (X_1/kappa)^n / (1 + (X_1/kappa)^n)]
% 				    R12:		 Z_1					--> 	phi				[delta_c*Z_1]
% 				    R13:		 Z_2					--> 	phi				[delta_c*Z_2]
% 				    R14:		 Z_3					--> 	phi				[delta_c*Z_3]
% 				    R15:		 Z_4					--> 	phi				[delta_c*Z_4]

%% Extract Parameters
k_1 = Parameters.k_1;
gamma_1 = Parameters.gamma_1;
gamma_2 = Parameters.gamma_2;
mu = Parameters.mu;
theta = Parameters.theta;
k = Parameters.k;
delta = Parameters.delta;
mu_0 = Parameters.mu_0;
alpha_1 = Parameters.alpha_1;
alpha_2 = Parameters.alpha_2;
eta = Parameters.eta;
eta_0 = Parameters.eta_0;
kappa = Parameters.kappa;
n = Parameters.n;
delta_c = Parameters.delta_c;

%% Extract State Variables
X_1 = x(1);
X_2 = x(2);
Z_1 = x(3);
Z_2 = x(4);
Z_3 = x(5);
Z_4 = x(6);

%% Propensities
a = [ ...
		gamma_1*X_1; ...
		gamma_2*X_2; ...
		k_1*X_1; ...
		mu; ...
		theta*X_2; ...
		eta*Z_1*Z_2; ...
		mu_0; ...
		alpha_1*Z_3 + alpha_2*X_2; ...
		eta_0*Z_3*Z_4; ...
		k*Z_1; ...
		(delta*X_2 + alpha_1*Z_3 + alpha_2*X_2) * (X_1/kappa)^n / (1 + (X_1/kappa)^n); ...
		delta_c*Z_1; ...
		delta_c*Z_2; ...
		delta_c*Z_3; ...
		delta_c*Z_4 ...
	];
end