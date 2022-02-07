%% Clear Workspace
close all;
clear;
clc;

addpath('Common Functions');

Load = 0;
if Load == 0
%% Fixed Parameters
RunTime = tic;
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

Parameters.k = 0.015;
Parameters.mu = 300 * 1e-2;
Parameters.theta = 0.3;
Parameters.eta = 1000 / 1e-2;
Parameters.kappa = 1 * 1e-6;
Parameters.kappa_0 = 1000 * 1e-2;
Parameters.alpha_0 = 100 * 1e-2;
Parameters.gamma_0 = 0.2;
Parameters.n = 1;
Parameters.n_0 = 1;
Parameters.delta_c = 0;

%% Setpoint, Disturbance and Settling Time Threshold
r = 1000 * 1e-2;
Factor_Disturbance = 2;
r_Disturbed = r * Factor_Disturbance;
Threshold = 0.01;
PerformanceWeights = [1, 1, 300];

%% Swept Parameters
N = 30;  % In paper: N = 300
delta_Sweep = linspace(0, 0.1, N);         
delta_0_Sweep = linspace(0, 4, N); 

%% Functions
StoichiometryMatrix = StoichiometryMatrix_APIDF3_Deg_Star();
PropensityFunction = @PropensityFunction_APIDF3_Deg_Star;
SupportingInput = @SupportingInput_Star;
JacobiansPlant = @Jacobians_Star;
FixedPoint = @APIDF3_Deg_FixedPoint;
Jacobians = @APIDF3_Deg_Jacobians;
Params2Gains = @APIDF3_Deg_Params2Gains;
OutputIndex = 6;

%% Simulation Settings
N_Simulation = 1000;
Solver = 'ODE15s';
t_1 = 3000;
t_f = 4000;

%% delta-delta_0 Sweep Simulations
PerformanceIndex = zeros(N,N);
K_I = zeros(N,N);
K_P = zeros(N,N);
K_D = zeros(N,N);
for i = 1 : N
    Parameters.delta_0 = delta_0_Sweep(i);
    for j = 1 : N
        Parameters.delta = delta_Sweep(j);
        [i, j]
    	IC = FixedPoint(Parameters, SupportingInput);
    	Parameters_Disturbed = Parameters;
        Parameters_Disturbed.mu = Parameters.mu * Factor_Disturbance;
        [A, ~, ~, ~] = Jacobians(Parameters_Disturbed, FixedPoint, SupportingInput, JacobiansPlant);
        lambda = eig(A);
        if real(lambda) <= 0
            [time_vector_1, X_1] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC, t_1, N_Simulation, Solver);
            [time_vector_2, X_2] = DSA(StoichiometryMatrix, PropensityFunction, Parameters_Disturbed, X_1(:,end), t_f-t_1, N_Simulation, Solver);
            PerformanceIndex(j,i) = Compute_PerformanceIndex(time_vector_2, X_2(OutputIndex,:), r_Disturbed, Threshold, PerformanceWeights);
            Gains = Params2Gains(Parameters_Disturbed, SupportingInput);
            K_I(j,i) = Gains.K_I;
            K_P(j,i) = Gains.K_P;
            K_D(j,i) = Gains.K_D;
     	else
            PerformanceIndex(j,i) = -10;
        end
    end
end
PerformanceIndex(PerformanceIndex == -10) = max(max(PerformanceIndex));
clear time_vector2 X_1 X_2 X

clear ans
RunTime = toc(RunTime);
save APIDF3_Deg_Star_PerformanceTracking_Sweep

else
load APIDF3_Deg_Star_PerformanceTracking_Sweep;
PerformanceIndex = PerformanceIndex/100;

%% K_P Contours
K_P_Contours = [-0.15, -0.1, -0.05, 0, 0.05];
delta_Contours = linspace(delta_Sweep(1), delta_Sweep(end), 100*N);
[u_bar, x_bar] = SupportingInput(Parameters_Disturbed, r_Disturbed);
delta_0_K_P_Contours = zeros(length(K_P_Contours), length(delta_Contours));
for i = 1 : length(K_P_Contours)
    T = (x_bar(1)/Parameters.kappa)^Parameters.n / (1 + (x_bar(1)/Parameters.kappa)^Parameters.n);
    sigma_6 = (Parameters.n_0 * Parameters.alpha_0 / r_Disturbed) * (r_Disturbed/Parameters.kappa_0)^Parameters.n_0 / (1 + (r_Disturbed/Parameters.kappa_0)^Parameters.n_0)^2;
    delta_0_K_P_Contours(i,:) = (Parameters.gamma_0/sigma_6) * (delta_Contours - K_P_Contours(i)/T);
    delta_0_K_P_Contours(i,delta_0_K_P_Contours(i,:) < 0) = NaN;
end

%% K_D Contours
K_D_Contours = [0.1, 0.2, 0.3, 0.4];
T = (x_bar(1)/Parameters.kappa)^Parameters.n / (1 + (x_bar(1)/Parameters.kappa)^Parameters.n);
delta_0_K_D_Contours = Parameters.gamma_0 * K_D_Contours / T;

%% Nominal Parameters: I, PI, ID, PID
delta_I_Nominal = 0;
delta_0_I_Nominal = 0;

delta_PID_Nominal_1 = 0.03;
delta_0_PID_Nominal_1 = 1.6;

% delta_PID_Nominal_2 = 0.0187;
% delta_0_PID_Nominal_2 = 2.1;
delta_PID_Nominal_2 = 0.05;
delta_0_PID_Nominal_2 = 2.8;

% PerformanceIndex(PerformanceIndex == min(min(PerformanceIndex))) = max(max(PerformanceIndex));
% [i, j] = find(PerformanceIndex == min(min(PerformanceIndex)), 1);
delta_PID_Nominal_3 = 0.01; % delta_Sweep(i);
delta_0_PID_Nominal_3 = 1.65; % delta_0_Sweep(j);

%% I Controller Nominal Simulation
Parameters.delta = delta_I_Nominal;
Parameters.delta_0 = delta_0_I_Nominal;
Parameters_Disturbed = Parameters;
IC = FixedPoint(Parameters, SupportingInput);
Parameters_Disturbed.mu = Parameters.mu * Factor_Disturbance;
[~, X_1] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC, t_1, N_Simulation, Solver);
[time_vector_2, X_2] = DSA(StoichiometryMatrix, PropensityFunction, Parameters_Disturbed, X_1(:,end), t_f-t_1, N_Simulation, Solver);
X = [X_1, X_2(:,2:end)];
x_AIF = X(OutputIndex,:);
SettlingTime_AIF_Nominal = Compute_SettlingTime(time_vector_2, X_2(OutputIndex,:), r_Disturbed, Threshold);
Gains_Nominal{1} = Params2Gains(Parameters_Disturbed, SupportingInput);
Gains_Nominal{1}
Parameters_I = Parameters_Disturbed;
clear time_vector2 X_1 X_2 X

%% PID Controller Nominal Simulation 1
Parameters.delta = delta_PID_Nominal_1;
Parameters.delta_0 = delta_0_PID_Nominal_1;
Parameters_Disturbed = Parameters;
IC = FixedPoint(Parameters, SupportingInput);
Parameters_Disturbed.mu = Parameters.mu * Factor_Disturbance;
[~, X_1] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC, t_1, N_Simulation, Solver);
[~, X_2] = DSA(StoichiometryMatrix, PropensityFunction, Parameters_Disturbed, X_1(:,end), t_f-t_1, N_Simulation, Solver);
X = [X_1, X_2(:,2:end)];
x_APIDF_1 = X(OutputIndex,:);
Gains_Nominal{2} = Params2Gains(Parameters_Disturbed, SupportingInput);
Gains_Nominal{2}
Parameters_PID_1 = Parameters_Disturbed;
clear time_vector2 X_1 X_2 X

%% PID Controller Nominal Simulation 2
Parameters.delta = delta_PID_Nominal_2;
Parameters.delta_0 = delta_0_PID_Nominal_2;
Parameters_Disturbed = Parameters;
IC = FixedPoint(Parameters, SupportingInput);
Parameters_Disturbed.mu = Parameters.mu * Factor_Disturbance;
[~, X_1] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC, t_1, N_Simulation, Solver);
[~, X_2] = DSA(StoichiometryMatrix, PropensityFunction, Parameters_Disturbed, X_1(:,end), t_f-t_1, N_Simulation, Solver);
X = [X_1, X_2(:,2:end)];
x_APIDF_2 = X(OutputIndex,:);
Gains_Nominal{3} = Params2Gains(Parameters_Disturbed, SupportingInput);
Gains_Nominal{3}
Parameters_PID_2 = Parameters_Disturbed;
clear time_vector2 X_1 X_2 X

%% PID Controller Nominal Simulation 3
Parameters.delta = delta_PID_Nominal_3;
Parameters.delta_0 = delta_0_PID_Nominal_3;
Parameters_Disturbed = Parameters;
IC = FixedPoint(Parameters, SupportingInput);
Parameters_Disturbed.mu = Parameters.mu * Factor_Disturbance;
[time_vector_1, X_1] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC, t_1, N_Simulation, Solver);
[time_vector_2, X_2] = DSA(StoichiometryMatrix, PropensityFunction, Parameters_Disturbed, X_1(:,end), t_f-t_1, N_Simulation, Solver);
X = [X_1, X_2(:,2:end)];
time_vector = [time_vector_1, t_1 + time_vector_2(2:end)];
x_APIDF_3 = X(OutputIndex,:);
RiseTime_APIDF_Nominal_3 = Compute_RiseTime(time_vector_2, X_2(OutputIndex,:), r_Disturbed, Threshold);
Gains_Nominal{4} = Params2Gains(Parameters_Disturbed, SupportingInput);
Gains_Nominal{4}
Parameters_PID_3 = Parameters_Disturbed;
clear time_vector2 X_1 X_2 X

%% Figure Settings
SS = 4;
Figure_Width = 8 * SS;
Figure_Height = 3.5 * SS;
FontSize = 5 * SS;
FontSize_Small = 3 * SS;
FontSize_Large = 6 * SS;
LineWidth = 0.65 * SS;
LineWidth_Thick = 1 * SS;
LineWidth_Thin = 0.5 * SS;
MarkerSize = 5 * SS;
IColor = [200, 200, 200]/255; % Blue
% IColor = [0, 130, 200]/255; % Blue
DColor1 = [148, 0, 211]/255; % Purple
DColor2 = [230, 25, 75]/255; % Red
DColor3 = [60, 180, 75]/255; % Green
TubeColor = [0.5, 0.5, 0.5]; % Gray
Opacity = 0.5;

%% Set Figure 1
Figure1_Name = 'APIDF3_Deg_Star';
Handle_Figure1 = figure();
    Handle_Figure1.Color = [1 1 1];
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [-30, 10 + 3.5*Figure_Height, 1.2*Figure_Width/2.5, 0.91*Figure_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];
    
%% Set Figure 4
Figure4_Name = 'APIDF3_Deg_NominalSimulations_Star';
Handle_Figure4 = figure();
    Handle_Figure4.Color = [1 1 1];
    Handle_Figure4.PaperUnits = 'centimeters';
    Handle_Figure4.Units = 'centimeters';
    Handle_Figure4.Position = [-30, 35.5, Figure_Width, 0.9*Figure_Height];
    Handle_Figure4.PaperPositionMode = 'auto';
    Handle_Figure4.PaperSize = [Handle_Figure4.PaperPosition(3), Handle_Figure4.PaperPosition(4)];

%% Axes for APIDF3 Settling Times and Contours
Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Position = [0.13, 0.08, 0.84, 0.88];
    Handle_Axis1.Box = 'on';
    Handle_Axis1.FontSize = FontSize;
    hold(Handle_Axis1, 'on');
    grid(Handle_Axis1, 'on');
    Handle_Axis1.XMinorGrid = 'off';
    Handle_Axis1.YMinorGrid = 'off';
    Handle_Axis1.Layer = 'top';
    Handle_Axis1.XScale = 'linear';
    Handle_Axis1.YScale = 'linear';
    Handle_Axis1.XLim = [delta_0_Sweep(1), delta_0_Sweep(end)];
    Handle_Axis1.YLim = [delta_Sweep(1), delta_Sweep(end)];
    Handle_Axis1.XLabel.String = '$\delta_0$';
    Handle_Axis1.XLabel.VerticalAlignment = 'top';
    Handle_Axis1.XLabel.Interpreter = 'latex';
    Handle_Axis1.YLabel.String = '$\delta$';
    Handle_Axis1.XLabel.Position(2) = 0;
    Handle_Axis1.YLabel.Position(1) = 0;
    Colormap = colormap(Handle_Axis1, 'turbo');
    colormap(Handle_Axis1, brighten(Colormap, 0));
    Handle_Axis1.YLabel.Interpreter = 'latex';
    Handle_Axis1.XTick = linspace(delta_0_Sweep(1), delta_0_Sweep(end), 6);
    Handle_Axis1.YTick = linspace(delta_Sweep(1), delta_Sweep(end), 6);
    Handle_Axis1.Position(3) = Handle_Axis1.Position(3) / 1.2;
    Handle_Axis1.Position(1) = Handle_Axis1.Position(1) / 1.2;
% Heat Map
Handle_IntensityMap1 = pcolor(Handle_Axis1, delta_0_Sweep, delta_Sweep, PerformanceIndex);
shading(Handle_Axis1, 'interp');
Handle_ColorBar = colorbar(Handle_Axis1);
    Handle_ColorBar.Position(1) = 0.85;
    Handle_ColorBar.Label.String = 'Performance Index';
% K_P Contours
plot(Handle_Axis1, delta_0_K_P_Contours, delta_Contours, 'Color', 'k', 'LineWidth', LineWidth_Thin/1.5, 'LineStyle', '-');
delta_0_Text = [3.8, 3.3, 2.4, 1.5, 0.4];
delta_Text = [0.08, 0.1, 0.1, 0.1, 0.1];
Rotation_Text = 68*[1, 1, 1, 1, 1];
Color_Text = ['k', 'w', 'k', 'w', 'w'];
for i = 1 : length(K_P_Contours) - 1
    Handle_Text = text(Handle_Axis1, delta_0_Text(i), delta_Text(i), ['$K_P = ', num2str(K_P_Contours(i)), '$']);
    Handle_Text.Color = Color_Text(i);
    Handle_Text.Interpreter = 'latex';
    Handle_Text.FontSize = FontSize * 0.9; 
    Handle_Text.VerticalAlignment = 'top';
    Handle_Text.HorizontalAlignment = 'right';
    Handle_Text.Rotation = Rotation_Text(i);
end
% K_D Contours
plot(Handle_Axis1,  Handle_Axis1.XLim, repmat(delta_0_K_D_Contours, 2, 1), 'Color', 'k', 'LineWidth', LineWidth_Thin/1.5, 'LineStyle', '-');
delta_0_Text = [0, 0, 0, 0];
delta_Text = K_D_Contours * Parameters.gamma_0;
Rotation_Text = 0*[1, 1, 1, 1];
Color_Text = ['w', 'w', 'w', 'w'];
for i = 1 : length(K_D_Contours)
    Handle_Text = text(Handle_Axis1, delta_0_Text(i), delta_Text(i), ['$K_D = ', num2str(K_D_Contours(i)), '$']);
    Handle_Text.Color = Color_Text(i);
    Handle_Text.Interpreter = 'latex';
    Handle_Text.FontSize = FontSize * 0.9; 
    Handle_Text.VerticalAlignment = 'bottom';
    Handle_Text.HorizontalAlignment = 'left';
    Handle_Text.Rotation = Rotation_Text(i);
end

% Nominal Values for Simulations
plot(Handle_Axis1, delta_0_I_Nominal, delta_I_Nominal, 'Marker', 'o', 'MarkerSize', MarkerSize/1.5, 'Color', 'k', 'MarkerFaceColor', IColor, 'LineWidth', LineWidth_Thin);
plot(Handle_Axis1, delta_0_PID_Nominal_1, delta_PID_Nominal_1, 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', 'k', 'MarkerFaceColor', DColor1, 'LineWidth', LineWidth_Thin);
plot(Handle_Axis1, delta_0_PID_Nominal_2, delta_PID_Nominal_2, 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', 'k', 'MarkerFaceColor', DColor2, 'LineWidth', LineWidth_Thin);
plot(Handle_Axis1, delta_0_PID_Nominal_3, delta_PID_Nominal_3, 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', 'k', 'MarkerFaceColor', DColor3, 'LineWidth', LineWidth_Thin);
    
%% Axis for Nominal Simulations
Handle_Axis4 = axes(Handle_Figure4);
    Handle_Axis4.Position = [0.09, 0.08, 0.88, 0.9];
    Handle_Axis4.Box = 'on';
    Handle_Axis4.FontSize = FontSize;
    hold(Handle_Axis4, 'on');
    grid(Handle_Axis4, 'on');
    Handle_Axis4.XMinorGrid = 'off';
    Handle_Axis4.YMinorGrid = 'off';
    Handle_Axis4.XLabel.String = 'Time';
    Handle_Axis4.YLabel.String = 'Output';
    Handle_Axis4.XScale = 'linear';
    Handle_Axis4.YScale = 'linear'; 
    Handle_Axis4.XLim = [0.99*t_1, t_f*3.8/4];
    Handle_Axis4.XTick = t_1 : 100 : t_f;
    Handle_Axis4.XTickLabel = 0 : 100 : t_f - t_1;
    Handle_Axis4.YLim = [0.9*r, 1.2*r_Disturbed];
    Handle_Axis4.XLabel.Position = [350 + t_1, Handle_Axis4.YLim(1)];
    Handle_Axis4.XLabel.VerticalAlignment = 'bottom';
% Rise Time
% plot(Handle_Axis4, t_1 + RiseTime_APIDF_Nominal_3 * [1, 1], [0, r_Disturbed * (1 + Threshold)], 'LineWidth', LineWidth_Thin, 'Color', DColor3, 'LineStyle', '--');
% Responses
Handle_R = plot(Handle_Axis4, [0, t_1, t_1, t_f], [r, r, r_Disturbed, r_Disturbed], 'Color', 'k', 'LineWidth', LineWidth_Thin, 'LineStyle', '--');
Handle_I = plot(Handle_Axis4, time_vector, x_AIF, 'LineWidth', LineWidth, 'Color', IColor);
Handle_PID_1 = plot(Handle_Axis4, time_vector, x_APIDF_1, 'LineWidth', LineWidth*1.5, 'Color', DColor1);
Handle_PID_3 = plot(Handle_Axis4, time_vector, x_APIDF_3, 'LineWidth', LineWidth*1.5, 'Color', DColor3);
Handle_PID_2 = plot(Handle_Axis4, time_vector, x_APIDF_2, 'LineWidth', LineWidth*1.5, 'Color', DColor2);
% Legends
Handle_Legend1 = legend(Handle_Axis4, [Handle_I, Handle_PID_1, Handle_PID_2, Handle_PID_3], {['I-Control, $K_I = ', num2str(round(Gains_Nominal{1}.K_I,3)), '$'], ...
                                       ['PID-Control, $K_I = ', num2str(round(Gains_Nominal{2}.K_I,3)), ',~K_P = ', num2str(round(Gains_Nominal{2}.K_P,3)), ', ~ K_D = ', num2str(round(Gains_Nominal{2}.K_D,3)), '$'], ...
                                       ['PID-Control, $K_I = ', num2str(round(Gains_Nominal{3}.K_I,3)), ',~K_P = ', num2str(round(Gains_Nominal{3}.K_P,3)), ', ~ K_D = ', num2str(round(Gains_Nominal{3}.K_D,3)), '$'], ...
                                       ['PID-Control, $K_I = ', num2str(round(Gains_Nominal{4}.K_I,3)), ',~K_P = ', num2str(round(Gains_Nominal{4}.K_P,3)), ', ~ K_D = ', num2str(round(Gains_Nominal{4}.K_D,3)), '$']});
    Handle_Legend1.FontSize = FontSize;
    Handle_Legend1.Interpreter = 'latex';
    Handle_Legend1.Location = 'Southeast';
 	Handle_Legend1.Position(2) = 2 * Handle_Legend1.Position(2);
    
%% Saving Figures
Save_Flag = 0;
if Save_Flag == 1  
    Handle_Figure1.Color = 'none';
    set(Handle_Figure1, 'InvertHardCopy', 'off');
    print(Handle_Figure1, Figure1_Name, '-dpdf', '-painters');
    
    Handle_Figure4.Color = 'none';
    set(Handle_Figure4, 'InvertHardCopy', 'off');
    print(Handle_Figure4, Figure4_Name, '-dpdf', '-painters');
end

end
