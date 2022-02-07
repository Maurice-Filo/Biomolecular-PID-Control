%% Clear Workspace
close all;
clear;
clc;
tic

Load = 1;
if Load == 0
%% Parameters that Control Accuracy and Computation Time
N = 300;                            % Discretization steps of Sweeping Proportional Control Parameters for Analytic Formula
N_Simulations = 5;                  % Discretization steps of Sweeping Proportional Control Parameters for Simulations (In Paper: N_Simulations = 100)
N_Trajectories = 400;               % Number of Trajectories used to compute statistics (transient) (In Paper: N_Trajectories = 10000)
tf = 100;                           % Final Time for Stochastic Transient Simulations
tf_Long = 1e3;                      % Final Time for Single Stochastic Simulation to reach Stationarity (In Paper: tf_Long = 2e4)
N_MaxRxnEvents = 2e5;               % Maximum Number of Events (used for N_Trajectory Simulation up to time = tf) (In Paper: N_MaxRxnEvents = 2e5)
N_MaxRxnEvents_Long = 1e6;          % Maximum Number of Events (used for 1 Trajectory Simulation up to time = tf_Long) (In Paper: N_MaxRxnEvents_Long = 4e7)

%% Fixed Parameters
RunTime = tic;
mu = 10;
theta = 2;
eta = 100;
k = 3;
kappa = 1;
kappa_1 = 1;
n = 1;
k_1 = 2;
gamma_1 = 2;
gamma_2 = 7;
r = mu/theta;

%% Swept Parameters
% Axis Limits
kappa_min = 0.2;
delta_max = 200;
Var_max = 13;
alpha_H = (gamma_1 * gamma_2 / k_1) * r * (1 + (r/kappa)^n);
alpha_max = 350;
alpha_vector = linspace(0, alpha_H, N);
kappa_vector = logspace(log10(kappa_min), 4, N);
delta_vector = linspace(0, delta_max, N);
alpha_vector_Simulations = linspace(0, alpha_H, N_Simulations);
kappa_vector_Simulations = logspace(log10(kappa_min), 4, N_Simulations);
delta_vector_Simulations = linspace(0, delta_max, N_Simulations);

%% Stochastic Simulation Settings
IC = zeros(4,1);
Portion = 1/5;                  % Portion of Trajectory considered to be transient (didn't reach stationarity)

%% Compute Approximate and Exact Stationary Variances
% AIF
clear Parameters
Parameters.mu = mu;
Parameters.theta = theta;
Parameters.k = k;
Parameters.gamma_1 = gamma_1;
Parameters.gamma_2 = gamma_2;
Parameters.k_1 = k_1;
Var_AIF = Approximate_AIF_Variance(Parameters);

% APIF Class 1 with Additive Inhibition
clear Parameters
Parameters.kappa = kappa;
Parameters.n = n;
Parameters.mu = mu;
Parameters.theta = theta;
Parameters.k = k;
Parameters.gamma_1 = gamma_1;
Parameters.gamma_2 = gamma_2;
Parameters.k_1 = k_1;
Parameters.eta = 1e4;
Parameters.beta = 0;
Parameters.gamma_c = 0;
Var_APIF1_Additive = zeros(length(alpha_vector), 1);
for i = 1 : length(alpha_vector)
    Parameters.alpha = alpha_vector(i);
    Var_APIF1_Additive(i) = Approximate_APIF1Additive_Variance(Parameters);
end
Var_APIF1_Additive_Simulations = zeros(length(alpha_vector_Simulations), 1);
for i = 1 : length(alpha_vector_Simulations)
    [1 i]
    Parameters.alpha = alpha_vector_Simulations(i);
    [T, X] = SSA(@PropensityFunction_APIDF_Additive_2LinearSpecies, StoichiometryMatrix_APIDF_Additive_2LinearSpecies(), Parameters, IC, 10*tf_Long, 10*N_MaxRxnEvents_Long);
    X_2 = X(2,:);
    N_0 = round(length(T) * Portion*2);
    Var_APIF1_Additive_Simulations(i) = (trapz(T(N_0:end), X_2(N_0:end).^2) / (T(end) - T(N_0))) - (trapz(T(N_0:end), X_2(N_0:end)) / (T(end) - T(N_0)))^2;
    clear T X X_2
end

% APIF Class 1 with Multiplicative Inhibition
clear Parameters
Parameters.n = n;
Parameters.mu = mu;
Parameters.theta = theta;
Parameters.k = k;
Parameters.gamma_1 = gamma_1;
Parameters.gamma_2 = gamma_2;
Parameters.k_1 = k_1;
Parameters.eta = 1e4;
Parameters.beta = 0;
Parameters.gamma_c = 0;
Var_APIF1_Multiplicative = zeros(length(kappa_vector), 1);
for i = 1 : length(kappa_vector)
    Parameters.kappa = kappa_vector(i);
    Var_APIF1_Multiplicative(i) = Approximate_APIF1Multiplicative_Variance(Parameters);
end
Var_APIF1_Multiplicative_Simulations = zeros(length(kappa_vector_Simulations), 1);
for i = 1 : length(kappa_vector_Simulations)
    [2 i]
    Parameters.kappa = kappa_vector_Simulations(i);
    [T, X] = SSA(@PropensityFunction_APIDF_Multiplicative_2LinearSpecies, StoichiometryMatrix_APIDF_Multiplicative_2LinearSpecies(), Parameters, IC, tf_Long, N_MaxRxnEvents_Long);
    X_2 = X(2,:);
    N_0 = round(length(T) * Portion);
    Var_APIF1_Multiplicative_Simulations(i) = (trapz(T(N_0:end), X_2(N_0:end).^2) / (T(end) - T(N_0))) - (trapz(T(N_0:end), X_2(N_0:end)) / (T(end) - T(N_0)))^2;
    clear T X X_2
end

% APIF Class 1 with Degradation Inhibition
clear Parameters
Parameters.n = n;
Parameters.mu = mu;
Parameters.theta = theta;
Parameters.k = k;
Parameters.gamma_1 = gamma_1;
Parameters.gamma_2 = gamma_2;
Parameters.k_1 = k_1;
Parameters.eta = 1e4;
Parameters.beta = 0;
Parameters.gamma_c = 0;
Parameters.kappa_1 = kappa_1;
Var_APIF1_Degradation = zeros(length(delta_vector), 1);
for i = 1 : length(delta_vector)
    Parameters.delta = delta_vector(i);
    Var_APIF1_Degradation(i) = Approximate_APIF1Degradation_Variance(Parameters);
end
Var_APIF1_Degradation_Simulations = zeros(length(delta_vector_Simulations), 1);
for i = 1 : length(delta_vector_Simulations)
    [3 i]
    Parameters.delta = delta_vector_Simulations(i);
    [T, X] = SSA(@PropensityFunction_APIDF_Degradation_2LinearSpecies, StoichiometryMatrix_APIDF_Degradation_2LinearSpecies(), Parameters, IC, tf_Long, N_MaxRxnEvents_Long);
    X_2 = X(2,:);
    N_0 = round(length(T) * Portion);
    Var_APIF1_Degradation_Simulations(i) = (trapz(T(N_0:end), X_2(N_0:end).^2) / (T(end) - T(N_0))) - (trapz(T(N_0:end), X_2(N_0:end)) / (T(end) - T(N_0)))^2;
    clear T X X_2
end

%% Unstable Additive
alpha_Unstable = 230;
clear Parameters
Parameters.alpha = alpha_Unstable;
Parameters.kappa = kappa;
Parameters.n = n;
Parameters.mu = mu;
Parameters.theta = theta;
Parameters.k = k;
Parameters.gamma_1 = gamma_1;
Parameters.gamma_2 = gamma_2;
Parameters.k_1 = k_1;
Parameters.eta = 1e4;
Parameters.beta = 0;
Parameters.gamma_c = 0;
T = cell(N_Trajectories,1);
X = cell(N_Trajectories,1);
parfor i = 1 : N_Trajectories
   [T{i}, X{i}] = SSA(@PropensityFunction_APIDF_Additive_2LinearSpecies, StoichiometryMatrix_APIDF_Additive_2LinearSpecies(), Parameters, IC, tf, N_MaxRxnEvents);
end
[time_vector, MeanX, VarX] = ComputeStatistics(T, X, 1000);
clear T X
0

%% Nominal Simulations
alpha_Nominal = 210;
kappa_Nominal = 1/3.2;
delta_Nominal = 150;

% Additive Inhibition
clear Parameters
Parameters.alpha = 0;
Parameters.kappa = kappa;
Parameters.n = n;
Parameters.mu = mu;
Parameters.theta = theta;
Parameters.k = k;
Parameters.gamma_1 = gamma_1;
Parameters.gamma_2 = gamma_2;
Parameters.k_1 = k_1;
Parameters.eta = 1e4;
Parameters.beta = 0;
Parameters.gamma_c = 0;
T = cell(N_Trajectories,1);
X = cell(N_Trajectories,1);
parfor i = 1 : N_Trajectories
   [T{i}, X{i}] = SSA(@PropensityFunction_APIDF_Additive_2LinearSpecies, StoichiometryMatrix_APIDF_Additive_2LinearSpecies(), Parameters, IC, tf, N_MaxRxnEvents);
end
[time_vector_AIF, ~, VarX_AIF] = ComputeStatistics(T, X, 1000);
clear T X
1

% Additive Inhibition
clear Parameters
Parameters.alpha = alpha_Nominal;
Parameters.kappa = kappa;
Parameters.n = n;
Parameters.mu = mu;
Parameters.theta = theta;
Parameters.k = k;
Parameters.gamma_1 = gamma_1;
Parameters.gamma_2 = gamma_2;
Parameters.k_1 = k_1;
Parameters.eta = 1e4;
Parameters.beta = 0;
Parameters.gamma_c = 0;
T = cell(N_Trajectories,1);
X = cell(N_Trajectories,1);
parfor i = 1 : N_Trajectories
   [T{i}, X{i}] = SSA(@PropensityFunction_APIDF_Additive_2LinearSpecies, StoichiometryMatrix_APIDF_Additive_2LinearSpecies(), Parameters, IC, tf, N_MaxRxnEvents);
end
[time_vector_Additive, ~, VarX_Additive] = ComputeStatistics(T, X, 1000);
clear T X
2

% Multiplicative Inhibition
clear Parameters
Parameters.kappa = kappa_Nominal;
Parameters.n = n;
Parameters.mu = mu;
Parameters.theta = theta;
Parameters.k = k;
Parameters.gamma_1 = gamma_1;
Parameters.gamma_2 = gamma_2;
Parameters.k_1 = k_1;
Parameters.eta = 1e4;
Parameters.beta = 0;
Parameters.gamma_c = 0;
T = cell(N_Trajectories,1);
X = cell(N_Trajectories,1);
parfor i = 1 : N_Trajectories
   [T{i}, X{i}] = SSA(@PropensityFunction_APIDF_Multiplicative_2LinearSpecies, StoichiometryMatrix_APIDF_Multiplicative_2LinearSpecies(), Parameters, IC, tf, N_MaxRxnEvents);
end
[time_vector_Multiplicative, ~, VarX_Multiplicative] = ComputeStatistics(T, X, 1000);
clear T X
3

% Degradation Inhibition
clear Parameters
Parameters.delta = delta_Nominal;
Parameters.n = n;
Parameters.mu = mu;
Parameters.theta = theta;
Parameters.k = k;
Parameters.gamma_1 = gamma_1;
Parameters.gamma_2 = gamma_2;
Parameters.k_1 = k_1;
Parameters.eta = 1e4;
Parameters.beta = 0;
Parameters.gamma_c = 0;
Parameters.kappa_1 = kappa_1;
T = cell(N_Trajectories,1);
X = cell(N_Trajectories,1);
parfor i = 1 : N_Trajectories
   [T{i}, X{i}] = SSA(@PropensityFunction_APIDF_Degradation_2LinearSpecies, StoichiometryMatrix_APIDF_Degradation_2LinearSpecies(), Parameters, IC, tf, N_MaxRxnEvents);
end
[time_vector_Degradation, ~, VarX_Degradation] = ComputeStatistics(T, X, 1000);
clear T X
4

RunTime = toc(RunTime);
save Results
else
load Results

% clear Parameters
% Parameters.n = n;
% Parameters.mu = mu;
% Parameters.theta = theta;
% Parameters.k = k;
% Parameters.gamma_1 = gamma_1;
% Parameters.gamma_2 = gamma_2;
% Parameters.k_1 = k_1;
% Parameters.eta = 1e4;
% Parameters.beta = 0;
% Parameters.gamma_c = 0;
% Parameters.kappa_1 = 1;
% Var_APIF1_Degradation = zeros(length(delta_vector), 1);
% for i = 1 : length(delta_vector)
%     Parameters.delta = delta_vector(i);
%     Var_APIF1_Degradation(i) = Approximate_APIF1Degradation_Variance(Parameters);
% end

%% Figure Settings
Colors = lines(10);
MyBlue = Colors(1,:); MyRed = Colors(2,:); MyGreen = Colors(5,:);
SS = 4; % Screen Scale
Figure_Width = 6 * SS;
Figure_Height = 3.5 * SS;
FontSize = 5 * SS;
FontSize_Small = 4 * SS;
LineWidth = 0.65 * SS;
LineWidth_Thick = 1 * SS;
MarkerSize = 3 * SS;

%% Set Figure 1
Figure1_Name = 'APIF1_Variance';
Handle_Figure1 = figure();
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [-30, 50, Figure_Width, Figure_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];
    
%% Set Figure 2
Figure2_Name = 'APIF1_NominalStochasticSimulations';
Handle_Figure2 = figure();
    Handle_Figure2.PaperUnits = 'centimeters';
    Handle_Figure2.Units = 'centimeters';
    Handle_Figure2.Position = [-30+Figure_Width, 50, Figure_Width, Figure_Height];
    Handle_Figure2.PaperPositionMode = 'auto';
    Handle_Figure2.PaperSize = [Handle_Figure2.PaperPosition(3), Handle_Figure2.PaperPosition(4)];

%% Axis for APIF1 (Additive) Stationary Variance
Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Position = [0.09, 0.18, 0.85, 0.73];
    Handle_Axis1.Box = 'on';
    Handle_Axis1.FontSize = FontSize;
    hold(Handle_Axis1, 'on');
    grid(Handle_Axis1, 'off');
    Handle_Axis1.YGrid = 'on';
    Handle_Axis1.XMinorGrid = 'off';
    Handle_Axis1.YMinorGrid = 'off';
    Handle_Axis1.XLabel.String = '$\alpha$';
    Handle_Axis1.YLabel.String = 'Var$_\pi[X_2]$';
    Handle_Axis1.XLabel.Interpreter = 'latex';
    Handle_Axis1.YLabel.Interpreter = 'latex';
    Handle_Axis1.XScale = 'linear';
    Handle_Axis1.YScale = 'linear'; 
    Handle_Axis1.XColor = MyBlue;
    Handle_Axis1.XLim = [0, alpha_max];
    Handle_Axis1.YLim = [0, Var_max];
    Handle_Axis1.XTick = [0, 100, 210, alpha_max];
    Handle_Axis1.XTickLabel = {'0', '100', '$\alpha_{TH} = 210$', num2str(alpha_max)};
    Handle_Axis1.TickLabelInterpreter = 'latex';
    Handle_Axis1.XLabel.Position(1:2) = [alpha_max/2, 0];
plot(Handle_Axis1, alpha_vector, Var_APIF1_Additive, 'LineWidth', LineWidth, 'Color', MyBlue, 'LineStyle', '--');
plot(Handle_Axis1, alpha_vector_Simulations, Var_APIF1_Additive_Simulations, 'LineWidth', LineWidth, 'Color', MyBlue, 'LineStyle', '-');
plot(Handle_Axis1, alpha_H*[1,1], Handle_Axis1.YLim, 'LineWidth', LineWidth/2, 'LineStyle', ':', 'Color', MyBlue);
plot(Handle_Axis1, alpha_Nominal, Var_APIF1_Additive_Simulations(find(alpha_vector_Simulations>=alpha_Nominal, 1)), 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', MyBlue, 'MarkerFaceColor', MyBlue);
plot(Handle_Axis1, 0, Var_APIF1_Additive_Simulations(1), 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', 'k', 'MarkerFaceColor', 'k');
Handle_Annotation1 = annotation(Handle_Figure1, 'textbox', 'LineStyle', 'none');
    Handle_Annotation1.String = 'Unstable';
    Handle_Annotation1.FontName = 'Helvetica';
    Handle_Annotation1.FontSize = FontSize/1.2;
    Handle_Annotation1.Position = [0.6 0.32 0 0];
    Handle_Annotation1.Color = MyBlue;
Handle_Annotation2 = annotation(Handle_Figure1, 'arrow', 'LineStyle', '-', 'Color', MyBlue);
    Handle_Annotation2.Position = [0.6, 0.26, 0.18, 0];
    
%% Axis for Unstable Additive
Handle_Axis1_Prime = axes(Handle_Figure1);
    Handle_Axis1_Prime.Position = [0.65, 0.6, 0.25, 0.25];
    Handle_Axis1_Prime.Box = 'on';
    Handle_Axis1_Prime.FontSize = FontSize;
    hold(Handle_Axis1_Prime, 'on');
    grid(Handle_Axis1_Prime, 'off');
    Handle_Axis1_Prime.XMinorGrid = 'off';
    Handle_Axis1_Prime.YMinorGrid = 'off';
    Handle_Axis1_Prime.XLabel.String = '$E[X_2]$';
    Handle_Axis1_Prime.YLabel.String = '$E[Z_2]$';
    Handle_Axis1_Prime.XLabel.Interpreter = 'latex';
    Handle_Axis1_Prime.YLabel.Interpreter = 'latex';
    Handle_Axis1_Prime.XTick = [];
    Handle_Axis1_Prime.YTick = [];
    Handle_Axis1_Prime.XScale = 'linear';
    Handle_Axis1_Prime.YScale = 'linear'; 
    Handle_Axis1_Prime.XLim = [0, 2.1*r];
    Handle_Axis1_Prime.XLabel.Position = [2*r,0];
    Handle_Axis1_Prime.YLabel.Position = [0,1.5];
    Handle_Axis1_Prime.Title.String = '$\alpha>\alpha_{TH}$';
    Handle_Axis1_Prime.Title.Interpreter = 'latex';
    Handle_Axis1_Prime.Title.Position = [r,1.95];
    Handle_Axis1_Prime.Title.Color = MyBlue;
    Handle_Axis1_Prime.XTick = r;
    Handle_Axis1_Prime.YLim = [0, 2];
    Handle_Axis1_Prime.XTickLabel = '$r$';
    Handle_Axis1_Prime.TickLabelInterpreter = 'latex';
plot(Handle_Axis1_Prime, MeanX(2,:), MeanX(4,:), 'LineWidth', LineWidth, 'Color', MyBlue, 'LineStyle', '-');
plot(Handle_Axis1_Prime, r*[1,1], [0,2], 'LineWidth', LineWidth/2, 'Color', 'k', 'LineStyle', '--');
% Handle_Patch = patch(Handle_Axis1_Prime, [MeanX(2,:) + VarX(2,:)/2, flip(MeanX(2,:) - VarX(2,:)/2)], [MeanX(4,:), flip(MeanX(4,:))], MyBlue);
%     Handle_Patch.EdgeColor = 'none';
%     Handle_Patch.FaceAlpha = 0.3;
    
%% Axis for APIF1 (Multiplicative) Stationary Variance
Handle_Axis2 = axes(Handle_Figure1);
    Handle_Axis2.Position = Handle_Axis1.Position;
    Handle_Axis2.Box = 'off';
    Handle_Axis2.FontSize = FontSize;
    hold(Handle_Axis2, 'on');
    grid(Handle_Axis2, 'off');
    Handle_Axis2.XMinorGrid = 'off';
    Handle_Axis2.YMinorGrid = 'off';
    Handle_Axis2.XLabel.String = '$\kappa^{-1}$';
    Handle_Axis2.XLabel.Interpreter = 'latex';
    Handle_Axis2.XScale = 'linear';
    Handle_Axis2.YScale = 'linear'; 
    Handle_Axis2.XAxisLocation = 'top';
    Handle_Axis2.XColor = MyRed;
    Handle_Axis2.Color = 'none';
    Handle_Axis2.XLim = [0, 1/kappa_min];
    Handle_Axis2.YLim = Handle_Axis1.YLim;
    Handle_Axis2.XTick = [0, 2, 4, 5];
    Handle_Axis2.YTick = [];
    Handle_Axis2.XLabel.Position(1:2) = [(1/kappa_min)/2, Handle_Axis1.YLim(2)];
plot(Handle_Axis2, 1./kappa_vector, Var_APIF1_Multiplicative, 'LineWidth', LineWidth, 'Color', MyRed, 'LineStyle', '--');
plot(Handle_Axis2, 1./kappa_vector_Simulations, Var_APIF1_Multiplicative_Simulations, 'LineWidth', LineWidth, 'Color', MyRed, 'LineStyle', '-');
plot(Handle_Axis2, 1/kappa_Nominal, Var_APIF1_Multiplicative_Simulations(find(kappa_vector_Simulations>=kappa_Nominal, 1)), 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', MyRed, 'MarkerFaceColor', MyRed);

%% Axis for APIF1 (Degradation) Stationary Variance
Handle_Axis3 = axes(Handle_Figure1);
    Handle_Axis3.Position = Handle_Axis1.Position;
    Handle_Axis3.Box = 'off';
    Handle_Axis2.FontSize = FontSize;
    hold(Handle_Axis3, 'on');
    grid(Handle_Axis3, 'off');
    Handle_Axis3.XMinorGrid = 'off';
    Handle_Axis3.YMinorGrid = 'off';
    Handle_Axis3.XScale = 'linear';
    Handle_Axis3.YScale = 'linear'; 
    Handle_Axis3.Visible = 'off';
    Handle_Axis3.XLim = [0, delta_max];
    Handle_Axis3.YLim = Handle_Axis1.YLim;
plot(Handle_Axis3, delta_vector, Var_APIF1_Degradation, 'LineWidth', LineWidth, 'Color', MyGreen, 'LineStyle', '--');
plot(Handle_Axis3, delta_vector_Simulations, Var_APIF1_Degradation_Simulations, 'LineWidth', LineWidth, 'Color', MyGreen, 'LineStyle', '-');
plot(Handle_Axis3, delta_Nominal, Var_APIF1_Degradation_Simulations(find(delta_vector_Simulations>=delta_Nominal, 1)), 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', MyGreen, 'MarkerFaceColor', MyGreen);

Handle_Axis3_Prime = axes(Handle_Figure1);
    Handle_Axis3_Prime.Position(1) = Handle_Axis3.Position(1);
    Handle_Axis3_Prime.Position(2) = 0.08;
    Handle_Axis3_Prime.Position(3) = Handle_Axis3.Position(3);
    Handle_Axis3_Prime.Position(4) = Handle_Axis3.Position(4);
    Handle_Axis3_Prime.FontSize = FontSize;
    Handle_Axis3_Prime.Color = 'none';
    Handle_Axis3_Prime.XColor = MyGreen;
    Handle_Axis3_Prime.YColor = 'none';
    Handle_Axis3_Prime.XTick = Handle_Axis3.XTick;
    Handle_Axis3_Prime.YTick = [];
    Handle_Axis3_Prime.XLabel.String = '$\delta$';
    Handle_Axis3_Prime.XLabel.Interpreter = 'latex';
    Handle_Axis3.YLabel.Interpreter = 'latex';
    Handle_Axis3_Prime.XLim = Handle_Axis3.XLim;
    Handle_Axis3_Prime.XTick = [0, 50, 150, 200];
    Handle_Axis3_Prime.XLabel.Position(1:2) = [delta_max/2, 0];
    
%% Legend
Handle_Axis_Legend1 = axes(Handle_Figure1);
    Handle_Axis_Legend1.Position = Handle_Axis1.Position;
    hold(Handle_Axis_Legend1, 'on');
patch(Handle_Axis_Legend1, -1, -1, MyBlue, 'EdgeColor' , 'none');
patch(Handle_Axis_Legend1, -1, -1, MyRed, 'EdgeColor' , 'none');
patch(Handle_Axis_Legend1, -1, -1, MyGreen, 'EdgeColor' , 'none');
Handle_Legend1 = legend(Handle_Axis_Legend1, {'Additive PI', 'Multiplicative PI', 'Degradation PI'});
    Handle_Legend1.FontSize = FontSize/1.2;
    Handle_Legend1.Position = [0.24, 0.82, 0, 0];
    Handle_Axis_Legend1.Visible = 'off';
    
Handle_Axis_Legend2 = axes(Handle_Figure1);
    Handle_Axis_Legend2.Position = Handle_Axis1.Position;
    hold(Handle_Axis_Legend2, 'on');
plot(Handle_Axis_Legend2, -1, -1, 'LineWidth', LineWidth, 'Color', [0.5, 0.5, 0.5]);
plot(Handle_Axis_Legend2, -1, -1, 'LineWidth', LineWidth, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--');
Handle_Legend2 = legend(Handle_Axis_Legend2, {'Exact', 'Approximate'});
    Handle_Legend2.FontSize = FontSize/1.2;
    Handle_Legend2.Position = [0.48, 0.82, 0, 0];
    Handle_Axis_Legend2.Visible = 'off';
    
%% Axis for Nominal Simulations
Handle_Axis4 = axes(Handle_Figure2);
    Handle_Axis4.Position = [0.09, 0.18, 0.85, 0.73];
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
    Handle_Axis4.XLim = [0, 80];
    Handle_Axis4.XLabel.Position(1:2) = [Handle_Axis4.XLim(2)/2, -1];
plot(Handle_Axis4, time_vector_AIF, VarX_AIF(2,:), 'LineWidth', LineWidth/2, 'Color', 'k');
plot(Handle_Axis4, time_vector_Additive, VarX_Additive(2,:), 'LineWidth', LineWidth, 'Color', MyBlue);
plot(Handle_Axis4, time_vector_Multiplicative, VarX_Multiplicative(2,:), 'LineWidth', LineWidth, 'Color', MyRed);
plot(Handle_Axis4, time_vector_Degradation, VarX_Degradation(2,:), 'LineWidth', LineWidth, 'Color', MyGreen);
Handle_Legend3 = legend(Handle_Axis4, {'Standalone I', 'Additive PI', 'Multiplicative PI', 'Degradation PI'});
    Handle_Legend3.FontSize = FontSize;
    
%% Saving Figures
Save_Flag = 0;
if Save_Flag == 1  
    Handle_Figure1.Color = 'none';
    set(Handle_Figure1, 'InvertHardCopy', 'off');
    print(Handle_Figure1, Figure1_Name,'-dpdf');
    
    Handle_Figure2.Color = 'none';
    set(Handle_Figure2, 'InvertHardCopy', 'off');
    print(Handle_Figure2, Figure2_Name,'-dpdf');
end
toc

end