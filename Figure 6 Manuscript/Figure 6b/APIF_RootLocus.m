%% Clear Workspace
close all;
clear;
clc;

Save_Flag = 0;

%% Plant Parameters
Parameters.gamma_1 = 0.1;
Parameters.gamma_2 = 0.1;
Parameters.k_1 = 0.1;

%% Fixed APIF Controller Parameters
delta_Nominal_1 = 35;
delta_Nominal_2 = 70;
delta_Nominal_3 = 150;
Parameters.mu = 10;
Parameters.theta = 2;
Parameters.eta = 1e4;
Parameters.kappa = 1e-5;
Parameters.k = 4;
Parameters.n = 1;
Parameters.delta_c = 0;

%% Setpoint, Disturbance and Settling Time Threshold
r = 5;
Factor_Disturbance = 2;
r_Disturbed = r * Factor_Disturbance;
Threshold = 0.02;

%% Swept Parameters
N = 10000;
delta_Sweep = linspace(30, 1500, N);   

%% Functions
StoichiometryMatrix = StoichiometryMatrix_APIF_Deg_GeneExp();
PropensityFunction = @PropensityFunction_APIF_Deg_GeneExp;
SupportingInput = @SupportingInput_GeneExp;
JacobiansPlant = @Jacobians_GeneExp;
FixedPoint = @APIF_Deg_FixedPoint;
Jacobians = @APIF_Deg_Jacobians;
Params2Gains = @APIF_Deg_Params2Gains;
OutputIndex = 2;

%% Compute Eigenvalues
lambda = zeros(4,N);
for i = 1 : N
    Parameters.delta = delta_Sweep(i);
    Parameters_Disturbed = Parameters;
    Parameters_Disturbed.mu = Parameters.mu * Factor_Disturbance;
    [A, ~, ~, ~] = Jacobians(Parameters_Disturbed, FixedPoint, SupportingInput, JacobiansPlant);
    lambda(:,i) = eig(A);
end

Parameters.delta = delta_Nominal_1;
[A, ~, ~, ~] = Jacobians(Parameters, FixedPoint, SupportingInput, JacobiansPlant);
lambda_Nominal_1 = eig(A);

Parameters.delta = delta_Nominal_2;
[A, ~, ~, ~] = Jacobians(Parameters, FixedPoint, SupportingInput, JacobiansPlant);
lambda_Nominal_2 = eig(A);

Parameters.delta = delta_Nominal_3;
[A, ~, ~, ~] = Jacobians(Parameters, FixedPoint, SupportingInput, JacobiansPlant);
lambda_Nominal_3 = eig(A);

%% Simulation Settings
N_Simulation = 1000;
Solver = 'ODE23s';
t_1 = 100;
t_f = 130;

%% Simulations
% AIF Controller 1
Parameters.delta = delta_Nominal_1;
Parameters_Disturbed = Parameters;
Parameters_Disturbed.mu = Parameters.mu * Factor_Disturbance;
IC = FixedPoint(Parameters, SupportingInput);
[time_vector_1, X_1] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC, t_1, N_Simulation, Solver);
[time_vector_2, X_2] = DSA(StoichiometryMatrix, PropensityFunction, Parameters_Disturbed, X_1(:,end), t_f-t_1, N_Simulation, Solver);
X = [X_1, X_2(:,2:end)];
time_vector = [time_vector_1, t_1 + time_vector_2(2:end)];
x_APIF_1 = X(OutputIndex, :);

% AIF Controller 2
Parameters.delta = delta_Nominal_2;
Parameters_Disturbed = Parameters;
Parameters_Disturbed.mu = Parameters.mu * Factor_Disturbance;
IC = FixedPoint(Parameters, SupportingInput);
[time_vector_1, X_1] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC, t_1, N_Simulation, Solver);
[time_vector_2, X_2] = DSA(StoichiometryMatrix, PropensityFunction, Parameters_Disturbed, X_1(:,end), t_f-t_1, N_Simulation, Solver);
X = [X_1, X_2(:,2:end)];
time_vector = [time_vector_1, t_1 + time_vector_2(2:end)];
x_APIF_2 = X(OutputIndex, :);

% AIF Controller 3
Parameters.delta = delta_Nominal_3;
Parameters_Disturbed = Parameters;
Parameters_Disturbed.mu = Parameters.mu * Factor_Disturbance;
IC = FixedPoint(Parameters, SupportingInput);
[time_vector_1, X_1] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC, t_1, N_Simulation, Solver);
[time_vector_2, X_2] = DSA(StoichiometryMatrix, PropensityFunction, Parameters_Disturbed, X_1(:,end), t_f-t_1, N_Simulation, Solver);
X = [X_1, X_2(:,2:end)];
time_vector = [time_vector_1, t_1 + time_vector_2(2:end)];
x_APIF_3 = X(OutputIndex, :);
        
%% Figure Settings
ScalingFactor = 1.25;
SS = 4;
Figure_Width = 8 * SS;
Figure_Height = 3 * SS;
FontSize = ScalingFactor*5 * SS;
FontSize_Small = ScalingFactor*3 * SS;
FontSize_Large = ScalingFactor*6 * SS;
LineWidth = ScalingFactor*0.65 * SS;
LineWidth_Thick = ScalingFactor*1 * SS;
LineWidth_Thin = ScalingFactor*0.25 * SS;
MarkerSize = ScalingFactor*7 * SS;
PIColor3 = [55,126,184]/255; % Blue
PIColor2 = [77,175,74]/255; % Green
PIColor1 = [228,26,28]/255; % Red

%% Set Figure 1
Figure1_Name = 'APIF_RootLocus_GeneExp';
Handle_Figure1 = figure();
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [-30, 10 + 3.5*Figure_Height, Figure_Width, Figure_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];

%% Axis for Root Locus (3D)
Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Position = [0.57, 0.12, 0.35, 0.84];
    Handle_Axis1.Box = 'on';
    Handle_Axis1.BoxStyle = 'full';
    Handle_Axis1.LineWidth = LineWidth_Thin;
    Handle_Axis1.FontSize = FontSize;
    hold(Handle_Axis1, 'on');
    grid(Handle_Axis1, 'on');
    Handle_Axis1.XMinorGrid = 'off';
    Handle_Axis1.YMinorGrid = 'off';
    Handle_Axis1.XLim = [0, 200];
    Handle_Axis1.YLim = [-0.3, 0.05];
    Handle_Axis1.ZLim = [-15, 15];
    Handle_Axis1.ZTick = [-15, 0, 15];
    Handle_Axis1.ZLabel.String = 'Imaginary Axis';
    Handle_Axis1.XLabel.Position(2) = Handle_Axis1.YLim(2);
    Handle_Axis1.ZLabel.Position(1) = Handle_Axis1.XLim(2);
    Handle_Axis1.ZLabel.Position(2) = Handle_Axis1.YLim(1)*1.25;
    Handle_Axis1.YTick = [-0.2, 0];
    Handle_Axis1.XTick = [0, 50, 100, 150, 200];
    Handle_Axis1.View = [103, 60];
  	set(Handle_Axis1, 'ColorScale', 'log');
    Colormap = colormap(Handle_Axis1, 'turbo');
Handle_Patch_Black = patch(Handle_Axis1, [Handle_Axis1.XLim, flip(Handle_Axis1.XLim)], -0.5 * (Parameters.gamma_1 + Parameters.gamma_2) * ones(1,4), [Handle_Axis1.ZLim(1) * [1, 1], Handle_Axis1.ZLim(2) * [1, 1]], 0.2*[1, 1, 1]);
    Handle_Patch_Black.FaceAlpha = 0.9;   
    Handle_Patch_Black.LineWidth = LineWidth_Thin;
    
Handle_Patch0 = patch(Handle_Axis1, [Handle_Axis1.XLim, flip(Handle_Axis1.XLim)], zeros(1,4), [Handle_Axis1.ZLim(1) * [1, 1], Handle_Axis1.ZLim(2) * [1, 1]], 'r');
    Handle_Patch0.FaceAlpha = 0.7;   
    Handle_Patch0.LineWidth = LineWidth_Thin;

Handle_Patch1 = patch(Handle_Axis1, delta_Nominal_1 * ones(1,4), [Handle_Axis1.YLim, flip(Handle_Axis1.YLim)], [Handle_Axis1.ZLim(1) * [1, 1], Handle_Axis1.ZLim(2) * [1, 1]], PIColor1);
    Handle_Patch1.FaceAlpha = 0.2;   
    Handle_Patch1.LineWidth = LineWidth_Thin;
    
Handle_Patch2 = patch(Handle_Axis1, delta_Nominal_2 * ones(1,4), [Handle_Axis1.YLim, flip(Handle_Axis1.YLim)], [Handle_Axis1.ZLim(1) * [1, 1], Handle_Axis1.ZLim(2) * [1, 1]], PIColor2);
    Handle_Patch2.FaceAlpha = 0.2;   
    Handle_Patch2.LineWidth = LineWidth_Thin;
    
Handle_Patch3 = patch(Handle_Axis1, delta_Nominal_3 * ones(1,4), [Handle_Axis1.YLim, flip(Handle_Axis1.YLim)], [Handle_Axis1.ZLim(1) * [1, 1], Handle_Axis1.ZLim(2) * [1, 1]], PIColor3);
    Handle_Patch3.FaceAlpha = 0.2;   
    Handle_Patch3.LineWidth = LineWidth_Thin;

for i = 1 : size(lambda,1)
	scatter3(Handle_Axis1, delta_Sweep, real(lambda(i,:)), imag(lambda(i,:)), MarkerSize, delta_Sweep, 'filled');
end
 
plot3(Handle_Axis1, delta_Nominal_1 * ones(length(lambda_Nominal_1), 1), real(lambda_Nominal_1), imag(lambda_Nominal_1), 'Marker', 'o', 'MarkerSize', MarkerSize/2, 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', PIColor1);
plot3(Handle_Axis1, delta_Nominal_2 * ones(length(lambda_Nominal_2), 1), real(lambda_Nominal_2), imag(lambda_Nominal_2), 'Marker', 'o', 'MarkerSize', MarkerSize/2, 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', PIColor2);
plot3(Handle_Axis1, delta_Nominal_3 * ones(length(lambda_Nominal_3), 1), real(lambda_Nominal_3), imag(lambda_Nominal_3), 'Marker', 'o', 'MarkerSize', MarkerSize/2, 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', PIColor3);

%% Axis for Root Locus (2D)
Handle_Axis3 = axes(Handle_Figure1);
    Handle_Axis3.Position = [0.07, 0.18, 0.38, 0.79];
    Handle_Axis3.Box = 'on';
    Handle_Axis3.BoxStyle = 'full';
    Handle_Axis3.LineWidth = LineWidth_Thin;
    Handle_Axis3.FontSize = FontSize;
    hold(Handle_Axis3, 'on');
    grid(Handle_Axis3, 'on');
    Handle_Axis3.XMinorGrid = 'off';
    Handle_Axis3.YMinorGrid = 'off';
    Handle_Axis3.XLim = [-0.3, 0.025];
    Handle_Axis3.YLim = [-15, 15];
    Handle_Axis3.XLabel.String = 'Real Axis';
    Handle_Axis3.XLabel.Position(2) = -18.5;
    Handle_Axis3.YLabel.String = 'Imaginary Axis';
    Handle_Axis3.YLabel.Position(1) = -0.334;
    Handle_Axis3.XTick = [-0.2, -(Parameters.gamma_1 + Parameters.gamma_2)/2, 0];
    Handle_Axis3.XTickLabel = {'-0.2', '', '0'};
    set(Handle_Axis3, 'ColorScale', 'log');
    colormap(Handle_Axis3, Colormap);
for i = 1 : size(lambda,1)
	scatter(Handle_Axis3, real(lambda(i,:)), imag(lambda(i,:)), MarkerSize, delta_Sweep, 'filled');
end
plot(Handle_Axis3, -0.5 * (Parameters.gamma_1 + Parameters.gamma_2) * ones(1,2), Handle_Axis3.YLim, 'Color', 'k', 'LineStyle', '-', 'LineWidth', LineWidth_Thin);
plot(Handle_Axis3, 0 * ones(1,2), Handle_Axis3.YLim, 'Color', 'r', 'LineStyle', '-', 'LineWidth', LineWidth_Thin);

plot(Handle_Axis3, real(lambda_Nominal_1), imag(lambda_Nominal_1), 'Marker', 'o', 'MarkerSize', MarkerSize/2, 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', PIColor1);
plot(Handle_Axis3, real(lambda_Nominal_2), imag(lambda_Nominal_2), 'Marker', 'o', 'MarkerSize', MarkerSize/2, 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', PIColor2);
plot(Handle_Axis3, real(lambda_Nominal_3), imag(lambda_Nominal_3), 'Marker', 'o', 'MarkerSize', MarkerSize/2, 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', PIColor3);

%% Annotations
Handle_Text2 = annotation(Handle_Figure1, 'textbox', [0.25, 0.07, 0.1, 0.1], 'String', '$-\frac{\gamma_1 + \gamma_2}{2}$', 'Interpreter', 'latex', 'FontSize', FontSize*1.2, 'EdgeColor', 'none'); 
Handle_Text1 = annotation(Handle_Figure1, 'textbox', [0.68, 0.06, 0.1, 0.1], 'String', '$-\frac{\gamma_1 + \gamma_2}{2}$', 'Interpreter', 'latex', 'FontSize', FontSize*1.2, 'EdgeColor', 'none'); 
Handle_YLabel = annotation(Handle_Figure1, 'textarrow', [0.74, 0.74], [0.03, 0.03], 'HeadStyle', 'none', 'LineStyle', 'none', 'String', 'Real Axis', 'FontSize', FontSize*1,  'TextRotation', -5); 
Handle_XLabel = annotation(Handle_Figure1, 'textarrow', [0.98, 0.98], [0.4, 0.4], 'HeadStyle', 'none', 'LineStyle', 'none', 'String', '$K_P = \delta$', 'FontSize', FontSize*1.2,  'TextRotation', 70, 'Interpreter', 'latex'); 

%% Set Figure 2
Figure2_Name = 'APIF_RootLocus_GeneExp_Simulation';
Handle_Figure2 = figure();
    Handle_Figure2.PaperUnits = 'centimeters';
    Handle_Figure2.Units = 'centimeters';
    Handle_Figure2.Position = [-30 + Figure_Width, 10 + 3.5*Figure_Height, Figure_Width/2 * 1.25, Figure_Height];
    Handle_Figure2.PaperPositionMode = 'auto';
    Handle_Figure2.PaperSize = [Handle_Figure2.PaperPosition(3), Handle_Figure2.PaperPosition(4)];

%% Axis for Simulations
Handle_Axis4 = axes(Handle_Figure2);
    Handle_Axis4.Position = [0.1, 0.18, 0.87, 0.74];
    Handle_Axis4.Box = 'on';
    Handle_Axis4.FontSize = FontSize;
    hold(Handle_Axis4, 'on');
    grid(Handle_Axis4, 'on');
    Handle_Axis4.XMinorGrid = 'off';
    Handle_Axis4.YMinorGrid = 'off';
    Handle_Axis4.XLim = [0.97*t_1, t_f];
    Handle_Axis4.XTick = t_1 : 10 : t_f;
    Handle_Axis4.XTickLabel = 0 : 10 : t_f - t_1;
    Handle_Axis4.XLabel.String = 'Time';
    Handle_Axis4.YLabel.String = 'Output';
    Handle_Axis4.YLim = [0.9*r, 1.1*r_Disturbed];
    Handle_Axis4.YLabel.Position(1) = 1.005*Handle_Axis4.YLabel.Position(1);
    Handle_Axis4.XLabel.Position(2) = 1.1*Handle_Axis4.XLabel.Position(2);
plot(Handle_Axis4, [0, t_1, t_1, t_f], [r, r, r_Disturbed, r_Disturbed], 'LineWidth', LineWidth_Thin, 'Color', 'k', 'LineStyle', '--');
plot(Handle_Axis4, time_vector, x_APIF_1, 'LineWidth', LineWidth, 'Color', PIColor1);
plot(Handle_Axis4, time_vector, x_APIF_2, 'LineWidth', LineWidth, 'Color', PIColor2);
plot(Handle_Axis4, time_vector, x_APIF_3, 'LineWidth', LineWidth, 'Color', PIColor3);
Handle_Legend = legend(Handle_Axis4, {'Set-Point', ['PI-Control $K_P = ', num2str(delta_Nominal_1), '$'], ['PI-Control $K_P = ', num2str(delta_Nominal_2), '$'], ['PI-Control $K_P = ', num2str(delta_Nominal_3), '$']});
Handle_Legend.Location = 'southeast';
Handle_Legend.Interpreter = 'latex';
Handle_Legend.FontSize = FontSize;

%% Save Figures
if Save_Flag == 1  
    Handle_Figure1.Color = 'none';
    set(Handle_Figure1, 'InvertHardCopy', 'off');
  	print(Handle_Figure1, Figure1_Name, '-dpdf', '-painters');
    Handle_Figure2.Color = 'none';
    set(Handle_Figure2, 'InvertHardCopy', 'off');
    print(Handle_Figure2, Figure2_Name, '-dpdf', '-painters');
end