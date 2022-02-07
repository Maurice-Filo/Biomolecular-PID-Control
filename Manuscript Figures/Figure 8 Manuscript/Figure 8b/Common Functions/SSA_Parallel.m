function [T_Stored, X_Stored, t_vector, MeanX, VarX, FinalDistributionX] = SSA_Parallel(Propensity, Stoichiometry, NetworkParameters, X0, tf, N_MaxRxnEvents, N_Total, N_Stored, N_Steps)
%% Uniform Time Grid
t_vector = linspace(0, tf, N_Steps);

%% Initialize Stored Trajectories
T_Stored = cell(N_Stored,1);
X_Stored = cell(N_Stored,1);
SumX = zeros(size(Stoichiometry,1),1);
SumX_Squared = zeros(size(Stoichiometry,1),1);
parfor i = 1 : N_Total
    % Simulate
    [T, X] = SSA(Propensity, Stoichiometry, NetworkParameters, X0, tf, N_MaxRxnEvents);
    % Store Trajectory
    if i <= N_Stored
        T_Stored{i} = T; X_Stored{i} = X;
    end
    % Resample at Time Grid
    X_Resampled = ResampleZOH(X, T, t_vector);
    % Compute the Sum
    SumX = SumX + X_Resampled;
    % Compute the Sum Squared
    SumX_Squared = SumX_Squared + X_Resampled.^2;
    % Store Final Value
    FinalDistributionX(:,i) = X_Resampled(:,end);
end
MeanX = SumX/N_Total;
VarX = SumX_Squared/(N_Total-1) - (N_Total/(N_Total-1)) * MeanX.^2;
end
  