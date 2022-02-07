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
delta_Nominal_1 = 300;
delta_Nominal_2 = 600;
Parameters.mu = 10;
Parameters.eta = 0.001;
Parameters.kappa = eps;
Parameters.n = 1;
Parameters.delta_c = 0;

%% Setpoint, Disturbance and Settling Time Threshold
r = 5;
Factor_Disturbance = 2;
r_Disturbed = r * Factor_Disturbance;
Threshold = 0.02;

%% Swept Parameters
N = 100;
delta_Sweep = linspace(0, 1000, N);         
beta_Sweep = linspace(0, 1000, N); 
theta_Sweep = beta_Sweep + Parameters.mu/r; 
k_Sweep = linspace(eps, 1000, N);

%% Functions
StoichiometryMatrix = StoichiometryMatrix_APIDF2_Deg_GeneExp();
PropensityFunction = @PropensityFunction_APIDF2_Deg_GeneExp;
SupportingInput = @SupportingInput_GeneExp;
JacobiansPlant = @Jacobians_GeneExp;
FixedPoint = @APIDF2_Deg_FixedPoint;
Jacobians = @APIDF2_Deg_Jacobians;
Params2Gains = @APIDF2_Deg_Params2Gains;
OutputIndex = 2;

%% Compute Eigenvalues
lambda = zeros(N,N,N,4);
lambda_Max = zeros(N,N,N);
for i = 1 : N
    i
    Parameters.beta = beta_Sweep(i);
    Parameters.theta = theta_Sweep(i);
    for j = 1 : N
        Parameters.delta = delta_Sweep(j);
        for l = 1 : N
            Parameters.k = k_Sweep(l);
            Parameters_Disturbed = Parameters;
            Parameters_Disturbed.mu = Parameters.mu * Factor_Disturbance;
            [A, ~, ~, ~] = Jacobians(Parameters_Disturbed, FixedPoint, SupportingInput, JacobiansPlant);
            EV = eig(A);
            lambda(i,j,l,:) = EV;
            EV_Max = EV(real(EV) == max(real(EV)));
            lambda_Max(i,j,l) = EV_Max(1);
        end
    end
end

% Parameters.delta = delta_Nominal_1;
% [A, ~, ~, ~] = Jacobians(Parameters, FixedPoint, SupportingInput, JacobiansPlant);
% lambda_Nominal_1 = eig(A);
% 
% Parameters.delta = delta_Nominal_2;
% [A, ~, ~, ~] = Jacobians(Parameters, FixedPoint, SupportingInput, JacobiansPlant);
% lambda_Nominal_2 = eig(A);

%% Simulation Settings
% N_Simulation = 1000;
% Solver = 'ODE23s';
% t_1 = 100;
% t_f = 150;

%% Simulations
% % AIF Controller 1
% Parameters.delta = delta_Nominal_1;
% Parameters_Disturbed = Parameters;
% Parameters_Disturbed.mu = Parameters.mu * Factor_Disturbance;
% IC = FixedPoint(Parameters, SupportingInput);
% [time_vector_1, X_1] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC, t_1, N_Simulation, Solver);
% [time_vector_2, X_2] = DSA(StoichiometryMatrix, PropensityFunction, Parameters_Disturbed, X_1(:,end), t_f-t_1, N_Simulation, Solver);
% X = [X_1, X_2(:,2:end)];
% time_vector = [time_vector_1, t_1 + time_vector_2(2:end)];
% x_APIF_1 = X(OutputIndex, :);
% 
% % AIF Controller 2
% Parameters.delta = delta_Nominal_2;
% Parameters_Disturbed = Parameters;
% Parameters_Disturbed.mu = Parameters.mu * Factor_Disturbance;
% IC = FixedPoint(Parameters, SupportingInput);
% [time_vector_1, X_1] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC, t_1, N_Simulation, Solver);
% [time_vector_2, X_2] = DSA(StoichiometryMatrix, PropensityFunction, Parameters_Disturbed, X_1(:,end), t_f-t_1, N_Simulation, Solver);
% X = [X_1, X_2(:,2:end)];
% time_vector = [time_vector_1, t_1 + time_vector_2(2:end)];
% x_APIF_2 = X(OutputIndex, :);
        
%% Figure Settings
SS = 4;
Figure_Width = 8 * SS;
Figure_Height = 3 * SS;
FontSize = 5 * SS;
FontSize_Small = 3 * SS;
FontSize_Large = 6 * SS;
LineWidth = 0.65 * SS;
LineWidth_Thick = 1 * SS;
LineWidth_Thin = 0.25 * SS;
MarkerSize = 7 * SS;
PIColor1 = [255, 127, 0]/255; % Orange
PIColor2 = [0, 130, 200]/255; % Blue
IColor = [0, 130, 200]/255; % Blue
DColor = [60, 180, 75]/255; % Green
TubeColor = [0.5, 0.5, 0.5]; % Gray
Opacity = 0.5;

%% Set Figure 1
Figure1_Name = 'APIDF4_RootLocus_GeneExp';
Handle_Figure1 = figure();
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [-30, 10 + 3.5*Figure_Height, Figure_Width/2.5, Figure_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];

%% Axis for Root Locus
Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Position = [0.1, 0.1, 0.8, 0.88];
    Handle_Axis1.Box = 'on';
    Handle_Axis1.BoxStyle = 'full';
    Handle_Axis1.LineWidth = LineWidth_Thin;
    Handle_Axis1.FontSize = FontSize;
    hold(Handle_Axis1, 'on');
    grid(Handle_Axis1, 'on');
    Handle_Axis1.XMinorGrid = 'off';
    Handle_Axis1.YMinorGrid = 'off';
    Handle_Axis1.XLim = [-10,1.5];
    Handle_Axis1.YLim = [-10, 10];
    Colormap = colormap(Handle_Axis1, [0, 0.6, 0; 1, 0, 0]);
scatter(Handle_Axis1, real(lambda(:)), imag(lambda(:)), MarkerSize, real(lambda(:)) >= 0, 'filled');
 
% plot3(Handle_Axis1, delta_Nominal_1 * ones(length(lambda_Nominal_1), 1), real(lambda_Nominal_1), imag(lambda_Nominal_1), 'Marker', 'o', 'MarkerSize', MarkerSize/2, 'Color', PIColor1, 'LineStyle', 'none', 'MarkerFaceColor', PIColor1);
% plot3(Handle_Axis1, delta_Nominal_2 * ones(length(lambda_Nominal_2), 1), real(lambda_Nominal_2), imag(lambda_Nominal_2), 'Marker', 'o', 'MarkerSize', MarkerSize/2, 'Color', PIColor2, 'LineStyle', 'none', 'MarkerFaceColor', PIColor2);
% 

%% Axis for Simulations
% Handle_Axis2 = axes(Handle_Figure1);
%     Handle_Axis2.Position = [0.07, 0.18, 0.4, 0.77];
%     Handle_Axis2.Box = 'on';
%     Handle_Axis2.FontSize = FontSize;
%     hold(Handle_Axis2, 'on');
%     grid(Handle_Axis2, 'on');
%     Handle_Axis2.XMinorGrid = 'off';
%     Handle_Axis2.YMinorGrid = 'off';
%     Handle_Axis2.XLim = [0.9*t_1, t_f];
%     Handle_Axis2.XLabel.String = 'Time';
%     Handle_Axis2.YLabel.String = 'Output';
%     Handle_Axis2.YLim = [0.9*r, 1.1*r_Disturbed];
% plot(Handle_Axis2, time_vector, x_APIF_1, 'LineWidth', LineWidth, 'Color', PIColor1);
% plot(Handle_Axis2, time_vector, x_APIF_2, 'LineWidth', LineWidth, 'Color', PIColor2);
% plot(Handle_Axis2, [0, t_1, t_1, t_f], [r, r, r_Disturbed, r_Disturbed], 'LineWidth', LineWidth_Thin, 'Color', 'k', 'LineStyle', '--');
% % Handle_Legend = legend(Handle_Axis2, 'Set-point', 'AIF Control');
% % Handle_Legend.Location = 'southeast';
% % Handle_Legend.FontSize = FontSize;
% 
% %% Annotations
% Handle_Line1 = annotation(Handle_Figure1, 'line', [0.674, 0.232], [0.1675, 0.738], 'LineWidth', LineWidth, 'LineStyle', ':', 'Color', PIColor1);
% Handle_Line2 = annotation(Handle_Figure1, 'line', [0.655, 0.371], [0.232, 0.77], 'LineWidth', LineWidth, 'LineStyle', ':', 'Color', PIColor2);
% 
% Handle_Text1 = annotation(Handle_Figure1, 'textbox', [0.85, 0.08, 0.1, 0.1], 'String', '$-\frac{\gamma_1 + \gamma_2}{2}$', 'Interpreter', 'latex', 'FontSize', FontSize*1.2, 'EdgeColor', 'none'); 
% Handle_YLabel = annotation(Handle_Figure1, 'textarrow', [0.88, 0.88], [0.06, 0.06], 'HeadStyle', 'none', 'LineStyle', 'none', 'String', 'Real', 'FontSize', FontSize*1.2,  'TextRotation', 10); 
% Handle_XLabel = annotation(Handle_Figure1, 'textarrow', [0.59, 0.59], [0.22, 0.22], 'HeadStyle', 'none', 'LineStyle', 'none', 'String', '$K_P = \delta$', 'FontSize', FontSize*1.2,  'TextRotation', -53, 'Interpreter', 'latex'); 

%% Save Figure
if Save_Flag == 1  
    Handle_Figure1.Color = 'none';
    set(Handle_Figure1, 'InvertHardCopy', 'off');
    print(Handle_Figure1, Figure1_Name, '-dpdf', '-painters');
end