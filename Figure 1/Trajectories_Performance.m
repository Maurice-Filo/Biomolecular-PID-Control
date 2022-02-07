%% Clear Workspace
close all;
clear;
clc;

Load  = 0;
if Load == 0
%% General Settings
RunTime = tic;
tf = 100;
t0 = 40;
x0 = zeros(4,1);
N = 1000;          
Solver = 'ODE23s';
StoichiometryMatrix = StoichiometryMatrix_APICType1Competitive_2SpeciesLinear();
PropensityFunction = @PropensityFunction_APICType1Competitive_2SpeciesLinear;
StoichiometryMatrix1 = StoichiometryMatrix_APICType1Additive_2SpeciesLinear();
PropensityFunction1 = @PropensityFunction_APICType1Additive_2SpeciesLinear;
r = 20;
N_SamplePaths = 100;     % In paper: N_SamplePaths = 1000
N_MaxRxnEvents = 5e4;
Factor_YLim = 2;

%% Ideal Response
NetworkParameters.mu = 100; 
NetworkParameters.theta = 5; 
NetworkParameters.eta = 100; 
NetworkParameters.delta = 0; 
NetworkParameters.k = 1; 
NetworkParameters.kappa = 0.15; 
NetworkParameters.n = 1;
NetworkParameters.k_1 = 2; 
NetworkParameters.gamma_1 = 5; 
NetworkParameters.gamma_2 = 0.1;
[time_vector, x] = DSA(StoichiometryMatrix, PropensityFunction, NetworkParameters, x0, tf, N, Solver);
y_Ideal = x(2,:);
clear x

%% Ideal Response Disturbed
NetworkParameters.mu = 100; 
NetworkParameters.theta = 5; 
NetworkParameters.eta = 100; 
NetworkParameters.delta = 0; 
NetworkParameters.k = 1; 
NetworkParameters.kappa = 0.15; 
NetworkParameters.n = 1;
NetworkParameters.k_1 = 2; 
NetworkParameters.gamma_1 = 5; 
NetworkParameters.gamma_2 = 0.1;
NetworkParameters_Disturbed = NetworkParameters;
NetworkParameters_Disturbed.gamma_2 = 1.6 * NetworkParameters.gamma_2;
[time_vector_Disturbance1, x] = DSA(StoichiometryMatrix, PropensityFunction, NetworkParameters, x0, t0, N, Solver);
y_Ideal_Disturbance = x(2,:);
[time_vector_Disturbance2, x] = DSA(StoichiometryMatrix, PropensityFunction, NetworkParameters_Disturbed, x(:,end), tf - t0, N, Solver);
y_Ideal_Disturbance = [y_Ideal_Disturbance, x(2,2:end)];
time_vector_Disturbance = [time_vector_Disturbance1, time_vector_Disturbance2(2:end) + t0];
clear NetworkParameters NetworkParameters_Disturbed
clear x time_vector_Disturbance1 time_vector_Disturbance2

%% Steady State Error
NetworkParameters.mu = 100; 
NetworkParameters.theta = 5; 
NetworkParameters.eta = 100; 
NetworkParameters.delta = 0; 
NetworkParameters.k = 1; 
NetworkParameters.kappa = 0.15; 
NetworkParameters.n = 1;
NetworkParameters.k_1 = 2; 
NetworkParameters.gamma_1 = 5; 
NetworkParameters.gamma_2 = 0.1;
NetworkParameters_Disturbed = NetworkParameters;
NetworkParameters_Disturbed.gamma_2 = 1.6 * NetworkParameters.gamma_2;
NetworkParameters_Disturbed.delta = 0.05;
[time_vector_Disturbance1, x] = DSA(StoichiometryMatrix, PropensityFunction, NetworkParameters, x0, t0, N, Solver);
y_SSError_Disturbance = x(2,:);
[time_vector_Disturbance2, x] = DSA(StoichiometryMatrix, PropensityFunction, NetworkParameters_Disturbed, x(:,end), tf - t0, N, Solver);
y_SSError_Disturbance = [y_SSError_Disturbance, x(2,2:end)];
time_vector_Disturbance_Error = [time_vector_Disturbance1, time_vector_Disturbance2(2:end) + t0];
clear NetworkParameters NetworkParameters_Disturbed
clear x time_vector_Disturbance1 time_vector_Disturbance2

%% Sustained Oscillations
NetworkParameters.mu = 100; 
NetworkParameters.theta = 5; 
NetworkParameters.eta = 100; 
NetworkParameters.delta = 0; 
NetworkParameters.k = 5; 
NetworkParameters.kappa = 3; 
NetworkParameters.n = 1;
NetworkParameters.k_1 = 2; 
NetworkParameters.gamma_1 = 5; 
NetworkParameters.gamma_2 = 0.1;
[~, x] = DSA(StoichiometryMatrix, PropensityFunction, NetworkParameters, x0, tf, N, Solver);
y_SustainedOscillations = x(2,:);
clear NetworkParameters
clear x

%% Divergent Response
NetworkParameters.mu = 0.8; 
NetworkParameters.theta = 5; 
NetworkParameters.eta = 100; 
NetworkParameters.delta = 0; 
NetworkParameters.k = 5; 
NetworkParameters.kappa = 100; 
NetworkParameters.n = 1;
NetworkParameters.k_1 = 2; 
NetworkParameters.gamma_1 = 0.005; 
NetworkParameters.gamma_2 = 0.005;
[~, x] = DSA(StoichiometryMatrix, PropensityFunction, NetworkParameters, x0, tf, N, Solver);
y_BlowUp = x(2,:);
clear NetworkParameters
clear x

%% Damped Oscillations
NetworkParameters.mu = 100; 
NetworkParameters.theta = 5; 
NetworkParameters.eta = 100; 
NetworkParameters.delta = 0; 
NetworkParameters.k = 5; 
NetworkParameters.kappa = 1.2; 
NetworkParameters.n = 1;
NetworkParameters.k_1 = 2; 
NetworkParameters.gamma_1 = 5; 
NetworkParameters.gamma_2 = 0.1;
[~, x] = DSA(StoichiometryMatrix, PropensityFunction, NetworkParameters, x0, tf, N, Solver);
y_DampedOscillations = x(2,:);
clear NetworkParameters
clear x

%% Slow Response
NetworkParameters.mu = 100; 
NetworkParameters.theta = 5; 
NetworkParameters.eta = 100; 
NetworkParameters.delta = 0; 
NetworkParameters.k = 5; 
NetworkParameters.kappa = 0.001; 
NetworkParameters.n = 1;
NetworkParameters.k_1 = 2; 
NetworkParameters.gamma_1 = 5; 
NetworkParameters.gamma_2 = 0.1;
[~, x] = DSA(StoichiometryMatrix, PropensityFunction, NetworkParameters, x0, tf, N, Solver);
y_Slow = x(2,:);
clear NetworkParameters
clear x

%% Overshoot
NetworkParameters.mu = 100; 
NetworkParameters.theta = 5; 
NetworkParameters.eta = 100; 
NetworkParameters.delta = 0; 
NetworkParameters.k = 2; 
NetworkParameters.kappa = 0.1; 
NetworkParameters.alpha = 1000;
NetworkParameters.n = 1;
NetworkParameters.k_1 = 2; 
NetworkParameters.gamma_1 = 5; 
NetworkParameters.gamma_2 = 0.1;
[~, x] = DSA(StoichiometryMatrix1, PropensityFunction1, NetworkParameters, x0, tf, N, Solver);
y_Overshoot = x(2,:);
clear NetworkParameters
clear x

%% Low Variance
NetworkParameters.mu = 100; 
NetworkParameters.theta = 5; 
NetworkParameters.eta = 100; 
NetworkParameters.delta = 0; 
NetworkParameters.k = 1; 
NetworkParameters.kappa = 0.05; 
NetworkParameters.n = 1;
NetworkParameters.k_1 = 2; 
NetworkParameters.gamma_1 = 5; 
NetworkParameters.gamma_2 = 0.1;
T_LowVariance = cell(N_SamplePaths, 1);
X_LowVariance = cell(N_SamplePaths, 1);
tic
for i = 1 : N_SamplePaths
    [T_LowVariance{i}, X_LowVariance{i}] = SSA(PropensityFunction, StoichiometryMatrix, NetworkParameters, x0, tf, N_MaxRxnEvents);
end
toc
[t_LowVariance, MeanX_LowVariance, VarX_LowVariance] = ComputeStatistics(T_LowVariance, X_LowVariance, N);
T_LowVariance = T_LowVariance{1};
X_LowVariance = X_LowVariance{1};

%% High Variance
NetworkParameters.mu = 100; 
NetworkParameters.theta = 5; 
NetworkParameters.eta = 100; 
NetworkParameters.delta = 0; 
NetworkParameters.k = 1; 
NetworkParameters.kappa = 5; 
NetworkParameters.n = 1;
NetworkParameters.k_1 = 2; 
NetworkParameters.gamma_1 = 5; 
NetworkParameters.gamma_2 = 0.1;
T_HighVariance = cell(N_SamplePaths, 1);
X_HighVariance = cell(N_SamplePaths, 1);
for i = 1 : N_SamplePaths
    [T_HighVariance{i}, X_HighVariance{i}] = SSA(PropensityFunction, StoichiometryMatrix, NetworkParameters, x0, tf, N_MaxRxnEvents);
end
[t_HighVariance, MeanX_HighVariance, VarX_HighVariance] = ComputeStatistics(T_HighVariance, X_HighVariance, N);
T_HighVariance = T_HighVariance{1};
X_HighVariance = X_HighVariance{1};

%% Save Results
RunTime = toc(RunTime);
save Simulations_Concept

else
load Simulations_Concept
%% Figure Settings
BackgroundColor = [230, 230, 230]/255;
Colors = lines(10);
MyBlue = Colors(1,:); MyRed = Colors(2,:); MyGreen = Colors(5,:);
SS = 4; % Screen Scale
Figure_Width = 5 * SS;
Figure_Height = 4 * SS;
FontSize = 8 * SS;
FontSize_Small = 4 * SS;
LineWidth_Thin = 0.4 * SS;
LineWidth = 0.8 * SS;
LineWidth_Thick = 1 * SS;
MarkerSize = 5 * SS;
% Figure 1
Figure1_Name = 'RPA Performance';
Handle_Figure1 = figure();
    Handle_Figure1.Color = [1, 1, 1];
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [5, 20, Figure_Width, Figure_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    set(Handle_Figure1, 'InvertHardCopy', 'off'); 
    
% Figure 2
Figure2_Name = 'Stability Performance';
Handle_Figure2 = figure();
    Handle_Figure2.Color = [1, 1, 1];
    Handle_Figure2.PaperUnits = 'centimeters';
    Handle_Figure2.Units = 'centimeters';
    Handle_Figure2.Position = [5, 20, Figure_Width, Figure_Height];
    Handle_Figure2.PaperPositionMode = 'auto';
    set(Handle_Figure2, 'InvertHardCopy', 'off');
    
% Figure 3
Figure3_Name = 'Transient Performance';
Handle_Figure3 = figure();
    Handle_Figure3.Color = [1, 1, 1];
    Handle_Figure3.PaperUnits = 'centimeters';
    Handle_Figure3.Units = 'centimeters';
    Handle_Figure3.Position = [5, 20, Figure_Width, Figure_Height];
    Handle_Figure3.PaperPositionMode = 'auto';
    set(Handle_Figure3, 'InvertHardCopy', 'off');
    
% Figure 4
Figure4_Name = 'Variance';
Handle_Figure4 = figure();
    Handle_Figure4.Color = [1, 1, 1];
    Handle_Figure4.PaperUnits = 'centimeters';
    Handle_Figure4.Units = 'centimeters';
    Handle_Figure4.Position = [5, 20, Figure_Width, Figure_Height];
    Handle_Figure4.PaperPositionMode = 'auto';
    set(Handle_Figure4, 'InvertHardCopy', 'off');
    
%% Plotting Figure 1
Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Box = 'on';
    Handle_Axis1.Position = [0.12 0.1 0.86 0.8];
    Handle_Axis1.FontSize = FontSize;
    hold(Handle_Axis1, 'on');
    grid(Handle_Axis1, 'off');
    Handle_Axis1.XMinorGrid = 'off';
    Handle_Axis1.YMinorGrid = 'off';
    Handle_Axis1.XLabel.String = '$t$';
    Handle_Axis1.YLabel.String = '$y(t)$';
    Handle_Axis1.Title.String = 'Robust Perfect Adaptation (RPA)';
    Handle_Axis1.Title.FontSize = FontSize;
    Handle_Axis1.YLim = [0, Factor_YLim*r];
    Handle_Axis1.XLim = [0, tf];
    Handle_Axis1.XTick = [];
    Handle_Axis1.YTick = r;
    Handle_Axis1.YTickLabel = '$r$';
    Handle_Axis1.XLabel.Interpreter = 'latex';
    Handle_Axis1.YLabel.Interpreter = 'latex';
    Handle_Axis1.TickLabelInterpreter = 'latex';
plot(Handle_Axis1, time_vector_Disturbance_Error, y_SSError_Disturbance, 'Color', Colors(2,:), 'LineWidth', LineWidth_Thick);
plot(Handle_Axis1, time_vector_Disturbance, y_Ideal_Disturbance, 'Color', Colors(1,:), 'LineWidth', LineWidth);
plot(Handle_Axis1, [time_vector(1), time_vector(N)], [r r], 'Color', 'k', 'LineWidth', LineWidth, 'LineStyle', '--');
Handle_Legend1 = legend(Handle_Axis1, {'RPA Lost', 'RPA Achieved', 'Set-point'});
Handle_Legend1.FontSize = FontSize*.8;
Handle_Legend1.Position(1) = 0.63;
Handle_Legend1.Position(2) = 0.7;

%% Plotting Figure 2
Handle_Axis2 = axes(Handle_Figure2);
    Handle_Axis2.Box = 'on';
    Handle_Axis2.Position = [0.12 0.1 0.86 0.8];
    Handle_Axis2.FontSize = FontSize;
    hold(Handle_Axis2, 'on');
    grid(Handle_Axis2, 'off');
    Handle_Axis2.XMinorGrid = 'off';
    Handle_Axis2.YMinorGrid = 'off';
    Handle_Axis2.XLabel.String = '$t$';
    Handle_Axis2.YLabel.String = '$y(t)$';
    Handle_Axis2.Title.String = 'Stability';
    Handle_Axis2.Title.Interpreter = 'tex';
    Handle_Axis2.YLim = [0, Factor_YLim*r];
    Handle_Axis2.XLim = [0, tf];
    Handle_Axis2.XTick = [];
    Handle_Axis2.YTick = r;
    Handle_Axis2.YTickLabel = '$r$';
    Handle_Axis2.XLabel.Interpreter = 'latex';
    Handle_Axis2.YLabel.Interpreter = 'latex';
    Handle_Axis2.TickLabelInterpreter = 'latex';
plot(Handle_Axis2, time_vector, y_BlowUp, 'Color', Colors(3,:), 'LineWidth', LineWidth_Thick);
plot(Handle_Axis2, time_vector, y_SustainedOscillations, 'Color', Colors(2,:), 'LineWidth', LineWidth_Thick);
plot(Handle_Axis2, time_vector, y_Ideal, 'Color', Colors(1,:), 'LineWidth', LineWidth);
plot(Handle_Axis2, [time_vector(1), time_vector(N)], [r r], 'Color', 'k', 'LineWidth', LineWidth, 'LineStyle', '--');
Handle_Legend2 = legend(Handle_Axis2, {'Divergent Response', 'Sustained Oscillations', 'Stable Fixed Point', 'Set-point'});
Handle_Legend2.FontSize = FontSize*.8;
Handle_Legend2.Position(1) = 0.47;
Handle_Legend2.Position(2) = 0.63;

%% Plotting Figure 3
Handle_Axis3 = axes(Handle_Figure3);
    Handle_Axis3.Box = 'on';
    Handle_Axis3.Position = [0.12 0.1 0.86 0.8];
    Handle_Axis3.FontSize = FontSize;
    hold(Handle_Axis3, 'on');
    grid(Handle_Axis3, 'off');
    Handle_Axis3.XMinorGrid = 'off';
    Handle_Axis3.YMinorGrid = 'off';
    Handle_Axis3.XLabel.String = '$t$';
    Handle_Axis3.YLabel.String = '$y(t)$';
    Handle_Axis3.Title.String = 'Transient Response';
    Handle_Axis3.Title.Interpreter = 'tex';
    Handle_Axis3.YLim = [0, Factor_YLim*r];
    Handle_Axis3.XLim = [0, tf];
    Handle_Axis3.XTick = [];
    Handle_Axis3.YTick = r;
    Handle_Axis3.YTickLabel = '$r$';
    Handle_Axis3.XLabel.Interpreter = 'latex';
    Handle_Axis3.YLabel.Interpreter = 'latex';
    Handle_Axis3.TickLabelInterpreter = 'latex';
plot(Handle_Axis3, time_vector, y_Slow, 'Color', Colors(3,:), 'LineWidth', LineWidth_Thick);
plot(Handle_Axis3, time_vector, y_DampedOscillations, 'Color', Colors(2,:), 'LineWidth', LineWidth_Thick);
plot(Handle_Axis3, time_vector, y_Overshoot, 'Color', Colors(4,:), 'LineWidth', LineWidth);
plot(Handle_Axis3, time_vector, y_Ideal, 'Color', Colors(1,:), 'LineWidth', LineWidth);
plot(Handle_Axis3, [time_vector(1), time_vector(N)], [r r], 'Color', 'k', 'LineWidth', LineWidth, 'LineStyle', '--');
Handle_Legend3 = legend(Handle_Axis3, {'Sluggish Response', 'Slowly Damped Oscillations', 'Overshoot', 'Smooth Transient', 'Set-point'});
Handle_Legend3.FontSize = FontSize*.8;
Handle_Legend3.Position(1) = 0.36;
Handle_Legend3.Position(2) = 0.56;

%% Plotting Figure 4
Handle_Axis4 = axes(Handle_Figure4);
    Handle_Axis4.Box = 'on';
    Handle_Axis4.Position = [0.12 0.1 0.86 0.8];
    Handle_Axis4.FontSize = FontSize;
    hold(Handle_Axis4, 'on');
    grid(Handle_Axis4, 'off');
    Handle_Axis4.XMinorGrid = 'off';
    Handle_Axis4.YMinorGrid = 'off';
    Handle_Axis4.XLabel.String = '$t$';
    Handle_Axis4.YLabel.String = '$Y(t)$';
    Handle_Axis4.Title.String = 'Stochastic Response';
    Handle_Axis4.Title.Interpreter = 'tex';
    Handle_Axis4.YLim = [0, Factor_YLim*r];
    Handle_Axis4.XLim = [0, tf];
    Handle_Axis4.XTick = [];
    Handle_Axis4.YTick = r;
    Handle_Axis4.YTickLabel = '$r$';
    Handle_Axis4.XLabel.Interpreter = 'latex';
    Handle_Axis4.YLabel.Interpreter = 'latex';
    Handle_Axis4.TickLabelInterpreter = 'latex';
Handle1 = patch(Handle_Axis4, [t_HighVariance, flip(t_HighVariance)], [MeanX_HighVariance(2,:) - 0.5*VarX_HighVariance(2,:), flip(MeanX_HighVariance(2,:) + 0.5*VarX_HighVariance(2,:))], Colors(2,:));
    Handle1.EdgeColor = 'none';
    Handle1.FaceAlpha = 0.3;
Handle2 = patch(Handle_Axis4, [t_LowVariance, flip(t_LowVariance)], [MeanX_LowVariance(2,:) - 0.5*VarX_LowVariance(2,:), flip(MeanX_LowVariance(2,:) + 0.5*VarX_LowVariance(2,:))], Colors(1,:));
    Handle2.EdgeColor = 'none';
    Handle2.FaceAlpha = 0.3;
plot(Handle_Axis4, t_HighVariance, MeanX_HighVariance(2,:), 'Color', Colors(2,:), 'LineWidth', LineWidth_Thick);
plot(Handle_Axis4, t_LowVariance, MeanX_LowVariance(2,:), 'Color', Colors(1,:), 'LineWidth', LineWidth_Thick);
plot(Handle_Axis4, T_HighVariance, X_HighVariance(2,:), 'Color', Colors(2,:), 'LineWidth', LineWidth_Thin);
plot(Handle_Axis4, T_LowVariance, X_LowVariance(2,:), 'Color', Colors(1,:), 'LineWidth', LineWidth_Thin);
Handle3 = plot(Handle_Axis4, [-1, -1], [-1, -1], 'Color', 0.5*[1,1,1], 'LineWidth', LineWidth_Thick);
Handle4 = plot(Handle_Axis4, [-1, -1], [-1, -1], 'Color', 0.5*[1,1,1], 'LineWidth', LineWidth_Thin);
Handle5 = plot(Handle_Axis4, [time_vector(1), time_vector(N)], [r r], 'Color', 'k', 'LineWidth', LineWidth, 'LineStyle', '--');
Handle_Legend4 = legend(Handle_Axis4, [Handle1, Handle2, Handle3, Handle4, Handle5], {'High Variance', 'Low Variance', 'Average', 'Example Trajectory', 'Set Point'});
Handle_Legend4.FontSize = FontSize*.6;
Handle_Legend4.Location = 'Southwest';
Handle_Legend4.Position(1) = 0.02;
Handle_Legend4.Position(2) = 0.02;

%% Saving Figures
Save_Flag = 0;
if Save_Flag == 1
    Handle_Figure1.Color = 'none';
    Handle_Figure2.Color = 'none';
    Handle_Figure3.Color = 'none';
    Handle_Figure4.Color = 'none';
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];
    Handle_Figure2.PaperUnits = 'centimeters';
    Handle_Figure2.PaperSize = [Handle_Figure2.PaperPosition(3), Handle_Figure2.PaperPosition(4)]; 
    Handle_Figure3.PaperUnits = 'centimeters';
    Handle_Figure3.PaperSize = [Handle_Figure3.PaperPosition(3), Handle_Figure3.PaperPosition(4)];
    Handle_Figure4.PaperUnits = 'centimeters';
    Handle_Figure4.PaperSize = [Handle_Figure4.PaperPosition(3), Handle_Figure4.PaperPosition(4)]; 
    print(Handle_Figure1, Figure1_Name,'-dpdf');
    print(Handle_Figure2, Figure2_Name,'-dpdf');
    print(Handle_Figure3, Figure3_Name,'-dpdf');
    print(Handle_Figure4, Figure4_Name,'-dpdf');
end

end

