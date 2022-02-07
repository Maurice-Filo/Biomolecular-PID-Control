%% Clear Workspace
close all;
clear;
clc;

Load = 0;
if Load == 0
%% Parameters that Control Accuracy and Computation Time
RunTime = tic;
N_Simulations = 1000;               % Discretization steps of Sweeping Proportional Control Parameters for Simulations
N_Trajectories = 64;                % Number of Trajectories used to compute statistics (transient)
tf = 200;                           % Final Time for Stochastic Transient Simulations
tf_Long = 2e4;                      % Final Time for Single Stochastic Simulation to reach Stationarity
N_MaxRxnEvents = 2.5e6;          	% Maximum Number of Events (used for N_Trajectory Simulation up to time = tf)
N_MaxRxnEvents_Long = 1e6;          % Maximum Number of Events (used for 1 Trajectory Simulation up to time = tf_Long) (In paper: N_MaxRxnEvents_Long = 1e8)
Portion = 1/5;                      % Portion of Trajectory considered to be transient (didn't reach stationarity)
N = 5;                              % Number of grid points for the swept parameter (In paper: N = 50)
Iterations = 1;               	    % Number of Outer Iterations for the batch stochastic Simulations: Each Iteration runs N_Trajectories using parfor (In paper: Iterations = 512)
Indeces_Nominal = round(linspace(1, N, 4));

%% Fixed Parameters
Parameters.gamma_1 = 2;
Parameters.gamma_2 = 7;
Parameters.k_1 = 2;

z_bar_3 = 1;
Parameters.k = 3;
Parameters.mu = 10;
Parameters.theta = 2;
Parameters.delta = 0;
Parameters.eta = 1e2;
Parameters.alpha_1 = 10;%3;
Parameters.mu_0 = 1000;
Parameters.eta_0 = 1e2;
Parameters.kappa = 1e-2;
Parameters.n = 1;
Parameters.delta_c = 0;

%% Setpoint
r = Parameters.mu/Parameters.theta;

%% Functions
StoichiometryMatrix = StoichiometryMatrix_APIDF4_Deg_GeneExp();
PropensityFunction = @PropensityFunction_APIDF4_Deg_GeneExp;
SupportingInput = @SupportingInput_GeneExp;
JacobiansPlant = @Jacobians_GeneExp;
FixedPoint = @APIDF4_Deg_FixedPoint;
Jacobians = @APIDF4_Deg_Jacobians;
Params2Gains = @APIDF4_Deg_Params2Gains;
Gains2Params = @APIDF4_Deg_Gains2Params;
OutputIndex = 2;

%% Swept K_D 
Parameters.alpha_2 = 0;
IC = zeros(6,1);
Gains = Params2Gains(Parameters, SupportingInput);
K_D_vector = linspace(0, 50, N);

%% K_D Sweep Simulations
StationaryVariance = zeros(N,1);
TStart0 = tic;
for i = 1 : N
    Gains.K_D = K_D_vector(i);
    [Parameters, Feasible] = Gains2Params(Gains, Parameters, r, SupportingInput);
    Parameters.mu_0 = Parameters.alpha_1*z_bar_3 + Parameters.alpha_2*r;
    [i, Feasible]
    [T, X] = SSA(PropensityFunction, StoichiometryMatrix(), Parameters, IC, tf_Long, N_MaxRxnEvents_Long);
	X_2 = X(OutputIndex,:);
  	N_0 = round(length(T) * Portion);
 	StationaryVariance(i) = (trapz(T(N_0:end), X_2(N_0:end).^2) / (T(end) - T(N_0))) - (trapz(T(N_0:end), X_2(N_0:end)) / (T(end) - T(N_0)))^2;
  	clear T X X_2
end
SimTime0 = toc(TStart0);

%% Nominal Derivative Gains for Transient Simulations
K_D_Nominal_1 = K_D_vector(Indeces_Nominal(1));
K_D_Nominal_2 = K_D_vector(Indeces_Nominal(2));
K_D_Nominal_3 = K_D_vector(Indeces_Nominal(3));
K_D_Nominal_4 = K_D_vector(Indeces_Nominal(4));
StationaryVariance_Nominal_1 = StationaryVariance(Indeces_Nominal(1));
StationaryVariance_Nominal_2 = StationaryVariance(Indeces_Nominal(2));
StationaryVariance_Nominal_3 = StationaryVariance(Indeces_Nominal(3));
StationaryVariance_Nominal_4 = StationaryVariance(Indeces_Nominal(4));

%% ID Controller Nominal Simulation 1
Gains.K_D = K_D_Nominal_1;
[Parameters_Nominal_1, Feasible_1] = Gains2Params(Gains, Parameters, r, SupportingInput);
Parameters_Nominal_1.mu_0 = Parameters_Nominal_1.alpha_1*z_bar_3 + Parameters_Nominal_1.alpha_2*r;
TStart1 = tic;
MeanY_Iteration = zeros(Iterations,N_Simulations-1);
VarY_Iteration = zeros(Iterations,N_Simulations-1);
FinalDistribution = zeros(Iterations,N_Trajectories);
for iter = 1 : Iterations
    [1, iter]
    T = cell(N_Trajectories,1);
    X = cell(N_Trajectories,1);
    Y = cell(N_Trajectories,1);
    parfor i = 1 : N_Trajectories
        [T{i}, X{i}] = SSA(PropensityFunction, StoichiometryMatrix, Parameters_Nominal_1, IC, tf, N_MaxRxnEvents);
        Y{i} = X{i}(OutputIndex,:);
        FinalDistribution(iter,i) = Y{i}(end);
    end
    [time_vector_1, MeanY_Iteration(iter,:), VarY_Iteration(iter,:)] = ComputeStatistics(T, Y, N_Simulations);
    clear T X Y
end
Mean_Nominal_1 = mean(MeanY_Iteration, 1);
Variance_Nominal_1 = mean(VarY_Iteration,1) + mean(MeanY_Iteration.^2,1) - Mean_Nominal_1.^2;
FinalDistribution1 = FinalDistribution(:);
SimTime1 = toc(TStart1);

%% ID Controller Nominal Simulation 2
Gains.K_D = K_D_Nominal_2;
[Parameters_Nominal_2, Feasible_2] = Gains2Params(Gains, Parameters, r, SupportingInput);
Parameters_Nominal_2.mu_0 = Parameters_Nominal_2.alpha_1*z_bar_3 + Parameters_Nominal_2.alpha_2*r;
TStart2 = tic;
MeanY_Iteration = zeros(Iterations,N_Simulations-1);
VarY_Iteration = zeros(Iterations,N_Simulations-1);
FinalDistribution = zeros(Iterations,N_Trajectories);
for iter = 1 : Iterations
    [2, iter]
    T = cell(N_Trajectories,1);
    X = cell(N_Trajectories,1);
    Y = cell(N_Trajectories,1);
    parfor i = 1 : N_Trajectories
        [T{i}, X{i}] = SSA(PropensityFunction, StoichiometryMatrix, Parameters_Nominal_2, IC, tf, N_MaxRxnEvents);
        Y{i} = X{i}(OutputIndex,:);
        FinalDistribution(iter,i) = Y{i}(end);
    end
    [time_vector_2, MeanY_Iteration(iter,:), VarY_Iteration(iter,:)] = ComputeStatistics(T, Y, N_Simulations);
    clear T X Y
end
Mean_Nominal_2 = mean(MeanY_Iteration, 1);
Variance_Nominal_2 = mean(VarY_Iteration,1) + mean(MeanY_Iteration.^2,1) - Mean_Nominal_2.^2;
FinalDistribution2 = FinalDistribution(:);
SimTime2 = toc(TStart2);

%% ID Controller Nominal Simulation 3
Gains.K_D = K_D_Nominal_3;
[Parameters_Nominal_3, Feasible_3] = Gains2Params(Gains, Parameters, r, SupportingInput);
Parameters_Nominal_3.mu_0 = Parameters_Nominal_3.alpha_1*z_bar_3 + Parameters_Nominal_3.alpha_2*r;
TStart3 = tic;
MeanY_Iteration = zeros(Iterations,N_Simulations-1);
VarY_Iteration = zeros(Iterations,N_Simulations-1);
FinalDistribution = zeros(Iterations,N_Trajectories);
for iter = 1 : Iterations
    [3, iter]
    T = cell(N_Trajectories,1);
    X = cell(N_Trajectories,1);
    Y = cell(N_Trajectories,1);
    parfor i = 1 : N_Trajectories
        [T{i}, X{i}] = SSA(PropensityFunction, StoichiometryMatrix, Parameters_Nominal_3, IC, tf, N_MaxRxnEvents);
        Y{i} = X{i}(OutputIndex,:);
        FinalDistribution(iter,i) = Y{i}(end);
    end
    [time_vector_3, MeanY_Iteration(iter,:), VarY_Iteration(iter,:)] = ComputeStatistics(T, Y, N_Simulations);
    clear T X Y
end
Mean_Nominal_3 = mean(MeanY_Iteration, 1);
Variance_Nominal_3 = mean(VarY_Iteration,1) + mean(MeanY_Iteration.^2,1) - Mean_Nominal_3.^2;
FinalDistribution3 = FinalDistribution(:);
SimTime3 = toc(TStart3);

%% ID Controller Nominal Simulation 4
Gains.K_D = K_D_Nominal_4;
[Parameters_Nominal_4, Feasible_4] = Gains2Params(Gains, Parameters, r, SupportingInput);
Parameters_Nominal_4.mu_0 = Parameters_Nominal_4.alpha_1*z_bar_3 + Parameters_Nominal_4.alpha_2*r;
TStart4 = tic;
MeanY_Iteration = zeros(Iterations,N_Simulations-1);
VarY_Iteration = zeros(Iterations,N_Simulations-1);
FinalDistribution = zeros(Iterations,N_Trajectories);
for iter = 1 : Iterations
    [4, iter]
    T = cell(N_Trajectories,1);
    X = cell(N_Trajectories,1);
    Y = cell(N_Trajectories,1);
    parfor i = 1 : N_Trajectories
        [T{i}, X{i}] = SSA(PropensityFunction, StoichiometryMatrix, Parameters_Nominal_4, IC, tf, N_MaxRxnEvents);
        Y{i} = X{i}(OutputIndex,:);
        FinalDistribution(iter,i) = Y{i}(end);
    end
    [time_vector_4, MeanY_Iteration(iter,:), VarY_Iteration(iter,:)] = ComputeStatistics(T, Y, N_Simulations);
    clear T X Y
end
Mean_Nominal_4 = mean(MeanY_Iteration, 1);
Variance_Nominal_4 = mean(VarY_Iteration,1) + mean(MeanY_Iteration.^2,1) - Mean_Nominal_4.^2;
FinalDistribution4 = FinalDistribution(:);
SimTime4 = toc(TStart4);

clear ans MeanY_Iteration VarY_Iteration iter
RunTime = toc(RunTime);
save APIDF4_Deg_GeneExp_Variance_Sweep

else
load APIDF4_Deg_GeneExp_Variance_Sweep;

%% Compute Biomolecular Parameters
N_Refined = 1000;
K_D_vector_Refined = linspace(K_D_vector(1), K_D_vector(end), N_Refined);
k_vector = zeros(N_Refined,1);
alpha_1_vector = zeros(N_Refined,1);
alpha_2_vector = zeros(N_Refined,1);
for i = 1 : N_Refined
    Gains.K_D = K_D_vector_Refined(i);
    [Parameters, Feasible] = Gains2Params(Gains, Parameters, r, SupportingInput);
    alpha_1_vector(i) = Parameters.alpha_1;
    alpha_2_vector(i) = Parameters.alpha_2;
end

%% Trajectories
N_VeryFewTrajectories = 16;

tic
T_1 = cell(N_VeryFewTrajectories,1);
Y_1 = cell(N_VeryFewTrajectories,1);
parfor i = 1 : N_VeryFewTrajectories
   [T_1{i}, X] = SSA(PropensityFunction, StoichiometryMatrix, Parameters_Nominal_1, IC, tf, N_MaxRxnEvents);
   Y_1{i} = X(OutputIndex,:);
end
toc

tic
T_3 = cell(N_VeryFewTrajectories,1);
Y_3 = cell(N_VeryFewTrajectories,1);
parfor i = 1 : N_VeryFewTrajectories
   [T_3{i}, X] = SSA(PropensityFunction, StoichiometryMatrix, Parameters_Nominal_3, IC, tf, N_MaxRxnEvents);
   Y_3{i} = X(OutputIndex,:);
end
toc

%% Figure Settings
SS = 4;
Figure_Width = 7.5 * SS;
Figure_Height = 3 * SS;
FontSize = 5 * SS;
FontSize_Small = 3 * SS;
FontSize_Large = 6 * SS;
LineWidth = 0.65 * SS;
LineWidth_Thick = 1 * SS;
LineWidth_Thin = 0.1 * SS;
MarkerSize = 5 * SS;
DColor1 = [148, 0, 211]/255; % Purple
DColor2 = [0, 130, 200]/255; % Blue
DColor3 = [60, 180, 75]/255; % Green
DColor4 = [230, 25, 75]/255; % Red
Opacity = 0.5;

%% Set Figure 1
Figure1_Name = 'APIDF4_Deg_GeneExp_BioParameters';
Handle_Figure1 = figure();
    Handle_Figure1.Color = [1 1 1];
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [-43, 10 + 3.5*Figure_Height, Figure_Width/2.5 * 1.1, Figure_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];
    
%% Set Figure 2
Figure2_Name = 'APIDF4_Deg_GeneExp_StationaryVariance';
Handle_Figure2 = figure();
    Handle_Figure2.Color = [1 1 1];
    Handle_Figure2.PaperUnits = 'centimeters';
    Handle_Figure2.Units = 'centimeters';
    Handle_Figure2.Position = [-43 + Figure_Width/2.5 * 1.1, 10 + 3.5*Figure_Height, Figure_Width/2.5 * 1.1, Figure_Height];
    Handle_Figure2.PaperPositionMode = 'auto';
    Handle_Figure2.PaperSize = [Handle_Figure2.PaperPosition(3), Handle_Figure2.PaperPosition(4)];
    
%% Set Figure 3
Figure3_Name = 'APIDF4_Deg_GeneExp_StochasticSimulations';
Handle_Figure3 = figure();
    Handle_Figure3.Color = [1 1 1];
    Handle_Figure3.PaperUnits = 'centimeters';
    Handle_Figure3.Units = 'centimeters';
    Handle_Figure3.Position = [-43 + 2*(Figure_Width/2.5 * 1.1), 10 + 3.5*Figure_Height, 2*(Figure_Width/2.5 * 1.1), Figure_Height];
    Handle_Figure3.PaperPositionMode = 'auto';
    Handle_Figure3.PaperSize = [Handle_Figure3.PaperPosition(3), Handle_Figure3.PaperPosition(4)];
    
%% Set Figure 4
Figure4_Name = 'APIDF4_Deg_NominalSimulations_GeneExp';
Handle_Figure4 = figure();
    Handle_Figure4.Color = [1 1 1];
    Handle_Figure4.PaperUnits = 'centimeters';
    Handle_Figure4.Units = 'centimeters';
    Handle_Figure4.Position = [-43 + 4*(Figure_Width/2.5 * 1.1), 10 + 3.5*Figure_Height, Figure_Width/2.5*1.1, Figure_Height];
    Handle_Figure4.PaperPositionMode = 'auto';
    Handle_Figure4.PaperSize = [Handle_Figure4.PaperPosition(3), Handle_Figure4.PaperPosition(4)];

%% Axes for BioParameters
Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Position = [0.18, 0.17, 0.6, 0.75];
    Handle_Axis1.Box = 'on';
    Handle_Axis1.BoxStyle = 'full';
    Handle_Axis1.LineWidth = LineWidth_Thin;
    Handle_Axis1.FontSize = FontSize;
    hold(Handle_Axis1, 'on');
    grid(Handle_Axis1, 'on');
    Handle_Axis1.XMinorGrid = 'off';
    Handle_Axis1.YMinorGrid = 'off';
    Handle_Axis1.XLabel.String = '$\alpha_2$';
    Handle_Axis1.YLabel.String = '$\alpha_1$';
    Handle_Axis1.XLabel.Interpreter = 'latex';
    Handle_Axis1.YLabel.Interpreter = 'latex';
    Handle_Axis1.XLim = [alpha_2_vector(1), alpha_2_vector(end)];
    Handle_Axis1.YLim = [0, 2*alpha_1_vector(end)];
    colormap(Handle_Axis1, 'turbo');
%         Handle_Axis1.Title.String = 'Biomolecular Parameters';
        Handle_ColorBar = colorbar(Handle_Axis1);
        Handle_ColorBar.Position(1) = 0.82;
        Handle_ColorBar.Label.String = '$K_D$';
        Handle_ColorBar.Label.Interpreter = 'latex';
scatter(Handle_Axis1, alpha_2_vector, alpha_1_vector, MarkerSize*3, linspace(K_D_vector_Refined(1), K_D_vector_Refined(N_Refined), N_Refined), 'filled');
plot(Handle_Axis1, Parameters_Nominal_1.alpha_2, Parameters_Nominal_1.alpha_1, 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', DColor1);
plot(Handle_Axis1, Parameters_Nominal_2.alpha_2, Parameters_Nominal_2.alpha_1, 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', DColor2);
plot(Handle_Axis1, Parameters_Nominal_3.alpha_2, Parameters_Nominal_3.alpha_1, 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', DColor3);
plot(Handle_Axis1, Parameters_Nominal_4.alpha_2, Parameters_Nominal_4.alpha_1, 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', DColor4);

%% Axes for Stationary Variance
Handle_Axis2 = axes(Handle_Figure2);
    Handle_Axis2.Position = [0.16, 0.16, 0.8, 0.76];
    Handle_Axis2.Box = 'on';
    Handle_Axis2.FontSize = FontSize;
    hold(Handle_Axis2, 'on');
    grid(Handle_Axis2, 'on');
    Handle_Axis2.Layer = 'top';
    Handle_Axis2.XMinorGrid = 'off';
    Handle_Axis2.YMinorGrid = 'off';
    Handle_Axis2.XScale = 'linear';
    Handle_Axis2.YScale = 'linear';
    Handle_Axis2.XLim = [K_D_vector(1), K_D_vector(end)];
    Handle_Axis2.XLabel.Interpreter = 'latex';
    Handle_Axis2.XLabel.String = '$K_D = \alpha_2/\alpha_1$';
    Handle_Axis2.YLabel.String = 'Var$_\pi[X_2]$';
    Handle_Axis2.YLabel.Interpreter = 'latex';
%     Handle_Axis2.Title.String = 'Stationary Variance';
plot(Handle_Axis2, K_D_vector, StationaryVariance, 'LineWidth', LineWidth, 'Color', 'k');
% Nominal Values for Simulations
plot(Handle_Axis2, K_D_Nominal_1, StationaryVariance_Nominal_1, 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', 'k', 'MarkerFaceColor', DColor1, 'LineWidth', LineWidth_Thin);
plot(Handle_Axis2, K_D_Nominal_2, StationaryVariance_Nominal_2, 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', 'k', 'MarkerFaceColor', DColor2, 'LineWidth', LineWidth_Thin);
plot(Handle_Axis2, K_D_Nominal_3, StationaryVariance_Nominal_3, 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', 'k', 'MarkerFaceColor', DColor3, 'LineWidth', LineWidth_Thin);
plot(Handle_Axis2, K_D_Nominal_4, StationaryVariance_Nominal_4, 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', 'k', 'MarkerFaceColor', DColor4, 'LineWidth', LineWidth_Thin);

%% Axes for Stochastic Trajectories
Handle_Axis31 = axes(Handle_Figure3);
    Handle_Axis31.Position = [0.08, 0.16, 0.49, 0.76];
    Handle_Axis31.Box = 'on';
    Handle_Axis31.FontSize = FontSize;
    hold(Handle_Axis31, 'on');
    grid(Handle_Axis31, 'on');
    Handle_Axis31.Layer = 'top';
    Handle_Axis31.XMinorGrid = 'off';
    Handle_Axis31.YMinorGrid = 'off';
    Handle_Axis31.XScale = 'linear';
    Handle_Axis31.YScale = 'linear';
    Handle_Axis31.XLabel.String = '$t$';
    Handle_Axis31.XLabel.Interpreter = 'latex';
    Handle_Axis31.YLabel.String = '$X_2(t)$';
    Handle_Axis31.YLabel.Interpreter = 'latex';
%     Handle_Axis31.Title.String = 'Stochastic Trajectories';
    Handle_Axis31.YLim = [0, 20];
    Handle_Axis31.XLim = [0, tf];
for i = 1 : N_VeryFewTrajectories
    plot(Handle_Axis31, T_1{i}, Y_1{i}, 'Color', [DColor1, 0.04], 'LineWidth', LineWidth_Thin);
    plot(Handle_Axis31, T_3{i}, Y_3{i}, 'Color', [DColor3, 0.05], 'LineWidth', LineWidth_Thin);
end
plot(Handle_Axis31, time_vector_1, Mean_Nominal_1, 'Color', DColor1, 'LineWidth', LineWidth_Thick);
plot(Handle_Axis31, time_vector_3, Mean_Nominal_3, 'Color', DColor3, 'LineWidth', LineWidth_Thick);

%% Axes for Final Distributions
Handle_Axis32 = axes(Handle_Figure3);
    Handle_Axis32.Position = [0.61, 0.16, 0.35, 0.76];
    Handle_Axis32.Box = 'on';
    Handle_Axis32.FontSize = FontSize;
    hold(Handle_Axis32, 'on');
    grid(Handle_Axis32, 'on');
    Handle_Axis32.Layer = 'top';
    Handle_Axis32.XMinorGrid = 'off';
    Handle_Axis32.YMinorGrid = 'off';
    Handle_Axis32.XScale = 'linear';
    Handle_Axis32.YScale = 'linear';
    Handle_Axis32.YLim = Handle_Axis31.YLim;
    Handle_Axis32.XLabel.String = 'Probability';
    Handle_Axis32.YTick = Handle_Axis31.YTick;
    Handle_Axis32.YTickLabel = [];
%     Handle_Axis32.Title.String = 'Stationary Distribution';
Handle_Histogram_1 = histogram(Handle_Axis32, FinalDistribution1, 'BinWidth', 1, 'FaceColor', DColor1, 'Orientation', 'Horizontal', 'Normalization', 'probability');
Handle_Histogram_3 = histogram(Handle_Axis32, FinalDistribution3, 'BinWidth', 1, 'FaceColor', DColor3, 'Orientation', 'Horizontal', 'Normalization', 'probability');
    Handle_Histogram_1.EdgeColor = 'none';
    Handle_Histogram_3.EdgeColor = 'none';
    
%% Axis for Nominal Simulations
Handle_Axis4 = axes(Handle_Figure4);
    Handle_Axis4.Position = [0.16, 0.15, 0.8, 0.76];
    Handle_Axis4.Box = 'on';
    Handle_Axis4.FontSize = FontSize;
    hold(Handle_Axis4, 'on');
    grid(Handle_Axis4, 'on');
    Handle_Axis4.XMinorGrid = 'off';
    Handle_Axis4.YMinorGrid = 'off';
    Handle_Axis4.XLabel.String = '$t$';
    Handle_Axis4.YLabel.String = 'Var$[X_2(t)]$';
    Handle_Axis4.XLabel.Interpreter = 'latex';
    Handle_Axis4.YLabel.Interpreter = 'latex';
    Handle_Axis4.XScale = 'linear';
    Handle_Axis4.YScale = 'linear'; 
    Handle_Axis4.XLim = [0, tf];
%     Handle_Axis4.Title.String = 'Variance Dynamics';
% Responses
Handle_ID_1 = plot(Handle_Axis4, time_vector_1, Variance_Nominal_1, 'LineWidth', LineWidth, 'Color', DColor1);
Handle_ID_2 = plot(Handle_Axis4, time_vector_2, Variance_Nominal_2, 'LineWidth', LineWidth, 'Color', DColor2);
Handle_ID_3 = plot(Handle_Axis4, time_vector_3, Variance_Nominal_3, 'LineWidth', LineWidth, 'Color', DColor3);
Handle_ID_4 = plot(Handle_Axis4, time_vector_4, Variance_Nominal_4, 'LineWidth', LineWidth, 'Color', DColor4);

%% Saving Figures
Save_Flag = 0;
if Save_Flag == 1  
    Handle_Figure1.Color = 'none';
    set(Handle_Figure1, 'InvertHardCopy', 'off');
  	print(Handle_Figure1, Figure1_Name, '-dpdf', '-painters');
    
    Handle_Figure2.Color = 'none';
    set(Handle_Figure2, 'InvertHardCopy', 'off');
  	print(Handle_Figure2, Figure2_Name, '-dpdf', '-painters');
    
    Handle_Figure3.Color = 'none';
    set(Handle_Figure3, 'InvertHardCopy', 'off');
  	print(Handle_Figure3, Figure3_Name, '-dpdf', '-painters');
    
    Handle_Figure4.Color = 'none';
    set(Handle_Figure4, 'InvertHardCopy', 'off');
  	print(Handle_Figure4, Figure4_Name, '-dpdf', '-painters');
end
end
