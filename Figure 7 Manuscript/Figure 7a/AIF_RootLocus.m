%% Clear Workspace
close all;
clear;
clc;

Save_Flag = 0;
addpath('Common Functions');

%% Plant Parameters
Parameters.gamma_F = 0.3;
Parameters.kappa_F = 1 * 1e-2;
Parameters.gamma_1 = 0.1;
Parameters.gamma_2 = 0.1;
Parameters.gamma_3 = 0.1;
Parameters.gamma_4 = 0.1;
Parameters.gamma_5 = 0.1;
Parameters.gamma_6 = 0.1;
Parameters.k_1 = 0.1;
Parameters.k_2 = 0.1;
Parameters.k_3 = 0.1;
Parameters.k_4 = 0.1;
Parameters.k_5 = 0.1;

%% Fixed AIF Controller Parameters
k_Nominal = 0.01;
Parameters.mu = 300 * 1e-2;
Parameters.theta = 0.3;
Parameters.eta = 1e2 / 1e-2;
Parameters.kappa = 1 * 1e-2;
Parameters.delta = 0;
Parameters.beta = 0;
Parameters.n = 1;
Parameters.delta_c = 0;
r = Parameters.mu / (Parameters.theta - Parameters.beta);

%% Swept Parameters
N = 50000;
k_Sweep = linspace(eps, 2, N);   

%% Functions
StoichiometryMatrix = StoichiometryMatrix_APIDF2_Deg_Star();
PropensityFunction = @PropensityFunction_APIDF2_Deg_Star;
SupportingInput = @SupportingInput_Star;
JacobiansPlant = @Jacobians_Star;
FixedPoint = @APIDF2_Deg_FixedPoint;
Jacobians = @APIDF2_Deg_Jacobians;
Params2Gains = @APIDF2_Deg_Params2Gains;
OutputIndex = 6;

%% Compute Eigenvalues
lambda = zeros(8,N);
for i = 1 : N
    Parameters.k = k_Sweep(i);
    [A, ~, ~, ~] = Jacobians(Parameters, FixedPoint, SupportingInput, JacobiansPlant);
    lambda(:,i) = eig(A);
end

Parameters.k = k_Nominal;
[A, ~, ~, ~] = Jacobians(Parameters, FixedPoint, SupportingInput, JacobiansPlant);
lambda_Nominal = eig(A);

%% Simulation Settings
N_Simulation = 1000;
Solver = 'ODE23s';
t_f = 1000;
IC = zeros(8,1);

%% Simulations
% Open Loop
Parameters.k = SupportingInput(Parameters, r);
[~, X] = DSA(StoichiometryMatrix_OL_Star(), @PropensityFunction_OL_Star, Parameters, IC(1:OutputIndex), t_f, N_Simulation, Solver);
x_OL = X(OutputIndex,:);
% AIF Controller
Parameters.k = k_Nominal;
[time_vector, X] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC, t_f, N_Simulation, Solver);
x_AIF = X(OutputIndex, :);

        
%% Figure Settings
SS = 4;
Figure_Width = 7 * SS;
Figure_Height = 3 * SS;
FontSize = 5 * SS;
FontSize_Small = 3 * SS;
FontSize_Large = 6 * SS;
LineWidth = 0.65 * SS;
LineWidth_Thick = 1 * SS;
LineWidth_Thin = 0.25 * SS;
MarkerSize = 7 * SS;
PColor = [230, 25, 75]/255; % Red
IColor = [0, 130, 200]/255; % Blue
DColor = [60, 180, 75]/255; % Green
HybridColor = DColor; % Green
TubeColor = [0.5, 0.5, 0.5]; % Gray
NeqColor = [0.5, 0.5, 0.5]; % Gray
PlantColor = [255, 127, 0]/255; % Orange
Opacity = 0.5;

%% Set Figure 1
Figure1_Name = 'AIF_RootLocus_Star';
Handle_Figure1 = figure();
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [-30 + Figure_Width/2, 10 + 3.5*Figure_Height, Figure_Width, Figure_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];

%% Set Figure 2
Figure2_Name = 'AIF_RootLocus_Star_2D';
Handle_Figure2 = figure();
    Handle_Figure2.PaperUnits = 'centimeters';
    Handle_Figure2.Units = 'centimeters';
    Handle_Figure2.Position = [-30, 10 + 3.5*Figure_Height, Figure_Width/2, Figure_Height];
    Handle_Figure2.PaperPositionMode = 'auto';
    Handle_Figure2.PaperSize = [Handle_Figure2.PaperPosition(3), Handle_Figure2.PaperPosition(4)];

%% Axis for Root Locus in 3D
Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Position = [0.09, 0.13, 0.35, 0.85];
    Handle_Axis1.Box = 'on';
    Handle_Axis1.BoxStyle = 'full';
    Handle_Axis1.LineWidth = LineWidth_Thin;
    Handle_Axis1.FontSize = FontSize;
    hold(Handle_Axis1, 'on');
    grid(Handle_Axis1, 'on');
    Handle_Axis1.XMinorGrid = 'off';
    Handle_Axis1.YMinorGrid = 'off';
    Handle_Axis1.XLim = [0, 0.05];
    Handle_Axis1.YLim = [-0.1, 0.05];
    Handle_Axis1.ZLim = [-0.1, 0.1];
    Handle_Axis1.XLabel.String = '$k$';
    Handle_Axis1.XLabel.Interpreter = 'latex';
    Handle_Axis1.YLabel.String = 'Real Axis';
    Handle_Axis1.YLabel.Rotation = -6;
    Handle_Axis1.ZLabel.String = 'Imaginary Axis';
    Handle_Axis1.XLabel.Position(2) = Handle_Axis1.YLim(2) * 1.8;
    Handle_Axis1.YLabel.Position(1) = 0.06;
    Handle_Axis1.YLabel.HorizontalAlignment = 'center';
    Handle_Axis1.ZLabel.Position(1) = 0.05;
    Handle_Axis1.ZLabel.Position(2) = Handle_Axis1.YLim(1)*1.3;
    Handle_Axis1.YTick = [-0.05, 0, 0.05];
    Handle_Axis1.XTick = [0, k_Nominal, 0.02, 0.03, 0.04];
    view(Handle_Axis1, [103, 60]);
    Colormap = colormap(Handle_Axis1, [0, 0.6, 0; 1, 0, 0]);
for i = 1 : size(lambda,1)
    scatter3(Handle_Axis1, k_Sweep, real(lambda(i,:)), imag(lambda(i,:)), MarkerSize, real(lambda(i,:)) >= 0, 'filled');
end
Handle_Patch1 = patch(Handle_Axis1, [Handle_Axis1.XLim, flip(Handle_Axis1.XLim)], zeros(1,4), [Handle_Axis1.ZLim(1) * [1, 1], Handle_Axis1.ZLim(2) * [1, 1]], 'r');
    Handle_Patch1.FaceAlpha = 0.1;   
    Handle_Patch1.LineWidth = LineWidth_Thin;

Handle_Patch2 = patch(Handle_Axis1, k_Nominal * ones(1,4), [Handle_Axis1.YLim, flip(Handle_Axis1.YLim)], [Handle_Axis1.ZLim(1) * [1, 1], Handle_Axis1.ZLim(2) * [1, 1]], IColor);
    Handle_Patch2.FaceAlpha = 0.1;   
    Handle_Patch2.LineWidth = LineWidth_Thin;
 
plot3(Handle_Axis1, k_Nominal * ones(length(lambda_Nominal), 1), real(lambda_Nominal), imag(lambda_Nominal), 'Marker', 'o', 'MarkerSize', MarkerSize/2, 'Color', IColor, 'LineStyle', 'none', 'MarkerFaceColor', IColor);

%% Axis for Simulations
Handle_Axis2 = axes(Handle_Figure1);
    Handle_Axis2.Position = [0.57, 0.18, 0.4, 0.77];
    Handle_Axis2.Box = 'on';
    Handle_Axis2.FontSize = FontSize;
    hold(Handle_Axis2, 'on');
    grid(Handle_Axis2, 'on');
    Handle_Axis2.XMinorGrid = 'off';
    Handle_Axis2.YMinorGrid = 'off';
    Handle_Axis2.XLim = [0, t_f];
    Handle_Axis2.XTick = linspace(0, t_f, 5);
    Handle_Axis2.XLabel.String = 'Time';
    Handle_Axis2.YLabel.String = 'Output';
    Handle_Axis2.YLim = [0, max([x_OL, x_AIF])];
plot(Handle_Axis2, Handle_Axis2.XLim, r * [1, 1], 'LineWidth', LineWidth_Thin, 'Color', 'k', 'LineStyle', '--');
plot(Handle_Axis2, time_vector, x_OL, 'LineWidth', LineWidth, 'Color', PlantColor);
plot(Handle_Axis2, time_vector, x_AIF, 'LineWidth', LineWidth, 'Color', IColor);
Handle_Legend = legend(Handle_Axis2, 'Set-point', 'Open Loop', 'I-Control, $K_I = 0.01$');
Handle_Legend.Interpreter = 'latex';
Handle_Legend.Location = 'southeast';
Handle_Legend.FontSize = FontSize;

%% Axis for Root Locus in 2D
Handle_Axis3 = axes(Handle_Figure2);
    Handle_Axis3.Position = [0.2, 0.16, 0.75, 0.8];
    Handle_Axis3.Box = 'on';
    Handle_Axis3.BoxStyle = 'full';
    Handle_Axis3.LineWidth = LineWidth_Thin;
    Handle_Axis3.FontSize = FontSize;
    hold(Handle_Axis3, 'on');
    grid(Handle_Axis3, 'on');
    Handle_Axis3.XMinorGrid = 'off';
    Handle_Axis3.YMinorGrid = 'off';
    Handle_Axis3.XLim = [-0.1, 0.05];
    Handle_Axis3.YLim = [-0.1, 0.1];
    Handle_Axis3.XLabel.String = 'Real Axis';
    Handle_Axis3.YLabel.String = 'Imaginary Axis';
    Handle_Axis3.XTick = [-0.05, 0, 0.05];
   	Colormap = colormap(Handle_Axis3, [0, 0.6, 0; 1, 0, 0]);
for i = 1 : size(lambda,1)
    scatter(Handle_Axis3, real(lambda(i,:)), imag(lambda(i,:)), MarkerSize, real(lambda(i,:)) >= 0, 'filled');
end
plot(Handle_Axis3, [0, 0], Handle_Axis3.YLim, 'LineWidth', LineWidth_Thin, 'Color', 'r');
plot(Handle_Axis3, real(lambda_Nominal), imag(lambda_Nominal), 'Marker', 'o', 'MarkerSize', MarkerSize/2, 'Color', IColor, 'LineStyle', 'none', 'MarkerFaceColor', IColor);


%% Save Figure
if Save_Flag == 1  
    Handle_Figure1.Color = 'none';
    set(Handle_Figure1, 'InvertHardCopy', 'off');
    print(Handle_Figure1, Figure1_Name, '-dpdf', '-painters');
    
    Handle_Figure2.Color = 'none';
    set(Handle_Figure2, 'InvertHardCopy', 'off');
    print(Handle_Figure2, Figure2_Name, '-dpdf', '-painters');
end