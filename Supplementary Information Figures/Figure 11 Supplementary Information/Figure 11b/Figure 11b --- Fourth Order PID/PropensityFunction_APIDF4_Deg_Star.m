function a = PropensityFunction_APIDF4_Deg_Star(x, Parameters)
% Propensity for Network APIDF4_Deg_Star
% 	 Species: 		 X = [X_1; X_2; X_3; X_4; X_5; X_6; Z_1; Z_2; Z_3; Z_4]
% 	 Reactions: 	R1:		 X_1					--> 	phi				[gamma_1*X_1]
% 				    R2:		 X_2					--> 	phi				[gamma_2*X_2]
% 				    R3:		 X_3					--> 	phi				[gamma_3*X_3]
% 				    R4:		 X_4					--> 	phi				[gamma_4*X_4]
% 				    R5:		 X_5					--> 	phi				[gamma_5*X_5]
% 				    R6:		 X_6					--> 	phi				[gamma_6*X_6]
% 				    R7:		 X_1					--> 	 X_1 +  X_2				[k_1*X_1]
% 				    R8:		 X_2					--> 	 X_2 +  X_3				[k_2*X_2]
% 				    R9:		 X_3					--> 	 X_3 +  X_4				[k_3*X_3]
% 				    R10:		 X_4					--> 	 X_4 +  X_5				[k_4*X_4]
% 				    R11:		 X_5					--> 	 X_5 +  X_6				[k_5*X_5]
% 				    R12:		 X_2					--> 	phi				[gamma_F*X_6*X_2 / (X_2 + kappa_F)]
% 				    R13:		phi					--> 	 Z_1				[mu]
% 				    R14:		phi					--> 	 Z_2				[theta*X_6]
% 				    R15:		 Z_1 +  Z_2					--> 	phi				[eta*Z_1*Z_2]
% 				    R16:		phi					--> 	 Z_3				[mu_0]
% 				    R17:		phi					--> 	 Z_4				[alpha_1*Z_3 + alpha_2*X_6]
% 				    R18:		 Z_3 +  Z_4					--> 	phi				[eta_0*Z_3*Z_4]
% 				    R19:		phi					--> 	 X_1				[k*Z_1]
% 				    R20:		 X_1					--> 	phi				[(delta*X_6 + alpha_1*Z_3 + alpha_2*X_6) * (X_1/kappa)^n / (1 + (X_1/kappa)^n)]
% 				    R21:		 Z_1					--> 	phi				[delta_c*Z_1]
% 				    R22:		 Z_2					--> 	phi				[delta_c*Z_2]
% 				    R23:		 Z_3					--> 	phi				[delta_c*Z_3]
% 				    R24:		 Z_4					--> 	phi				[delta_c*Z_4]

%% Extract Parameters
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
X_3 = x(3);
X_4 = x(4);
X_5 = x(5);
X_6 = x(6);
Z_1 = x(7);
Z_2 = x(8);
Z_3 = x(9);
Z_4 = x(10);

%% Propensities
a = [ ...
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
		gamma_F*X_6*X_2 / (X_2 + kappa_F); ...
		mu; ...
		theta*X_6; ...
		eta*Z_1*Z_2; ...
		mu_0; ...
		alpha_1*Z_3 + alpha_2*X_6; ...
		eta_0*Z_3*Z_4; ...
		k*Z_1; ...
		(delta*X_6 + alpha_1*Z_3 + alpha_2*X_6) * (X_1/kappa)^n / (1 + (X_1/kappa)^n); ...
		delta_c*Z_1; ...
		delta_c*Z_2; ...
		delta_c*Z_3; ...
		delta_c*Z_4 ...
	];
end