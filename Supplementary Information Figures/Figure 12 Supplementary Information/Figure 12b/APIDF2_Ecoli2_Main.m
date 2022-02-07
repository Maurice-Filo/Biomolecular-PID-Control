%% Clear Workspace
close all
clear
clc
Save_Flag = 0;
addpath('Common Functions');

%% Functions
StoichiometryMatrix = StoichiometryMatrix_APIDF2_Ecoli2();
PropensityFunction = @PropensityFunction_APIDF2_Ecoli2;
IdealSupportingInput = @IdealSupportingInput_Ecoli2;
Params2Gains = @APIDF2_Params2Gains;
OutputIndex = 2;

%% Plant Parameters
Parameters.gamma_1 = 0.028;             % [1/min]
Parameters.gamma_2 = 0.028;             % [1/min]
Parameters.k_1 = 1;                     % [1/min]

%% Controller Parameters
Parameters.delta_c = 0.028;             % [1/min]
Parameters.mu = 100;                   % [nM/min]
Parameters.theta = 1;                   % [1/min]
Parameters.eta = 0.05;                % [1/nM/min]
Parameters.k = 0.05;                       % [1/min]

Parameters.kappa = 10;                   % [nM]    
Parameters.delta = 0;           % [1/min]
Parameters.n = 1;                       % [Dimensionless]

Parameters.beta = 0;                    % [1/min]
r = Parameters.mu/(Parameters.theta - Parameters.beta);

%% Disturbance Parameters
gamma_2 = 0.028;

%% PID Parameters
% I-Control
mu_I = 100;
delta_I = 0;
beta_I = 0;
% PI-Control
mu_PI = 100;
delta_PI = 5 * 0.028;
beta_PI = 0;
% PID-Control
mu_PID = 50;
delta_PID = 5 * 0.028;
beta_PID = 0.5;

%% Simulation Settings
N_Simulation = 10000;
Solver = 'ODE15s';
t_1 = 8 * 60;
t_2 = 16 * 60;
t_f = 24 * 60;
TrackingFactor = 2;
DisturbanceFactor = 4;

%% Simulations
IC_1 = zeros(4,1);
% I-Control
Parameters.gamma_2 = gamma_2;
Parameters.mu = mu_I/TrackingFactor; Parameters.delta = delta_I; Parameters.beta = beta_I;
[time_I_1, X_I_1] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC_1, t_1, N_Simulation, Solver);
IC_2 = X_I_1(:,end);
Parameters.mu = mu_I;
[time_I_2, X_I_2] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC_2, t_2-t_1, N_Simulation, Solver);
IC_3 = X_I_2(:,end);
Parameters.gamma_2 = gamma_2*DisturbanceFactor;
[time_I_3, X_I_3] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC_3, t_f-t_2, N_Simulation, Solver);
time_I = [time_I_1, time_I_2(2:end) + t_1, time_I_3(2:end) + t_2];
X_I = [X_I_1, X_I_2(:,2:end), X_I_3(:,2:end)];
y_I = X_I(OutputIndex,:);
% PI-Control
Parameters.gamma_2 = gamma_2;
Parameters.mu = mu_PI/TrackingFactor; Parameters.delta = delta_PI; Parameters.beta = beta_PI;
[time_PI_1, X_PI_1] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC_1, t_1, N_Simulation, Solver);
IC_2 = X_PI_1(:,end);
Parameters.mu = mu_PI;
[time_PI_2, X_PI_2] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC_2, t_2-t_1, N_Simulation, Solver);
IC_3 = X_PI_2(:,end);
Parameters.gamma_2 = gamma_2*DisturbanceFactor;
[time_PI_3, X_PI_3] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC_3, t_f-t_2, N_Simulation, Solver);
time_PI = [time_PI_1, time_PI_2(2:end) + t_1, time_PI_3(2:end) + t_2];
X_PI = [X_PI_1, X_PI_2(:,2:end), X_PI_3(:,2:end)];
y_PI = X_PI(OutputIndex,:);
% PID-Control
Parameters.gamma_2 = gamma_2;
Parameters.mu = mu_PID/TrackingFactor; Parameters.delta = delta_PID; Parameters.beta = beta_PID;
[time_PID_1, X_PID_1] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC_1, t_1, N_Simulation, Solver);
IC_2 = X_PID_1(:,end);
Parameters.mu = mu_PID;
[time_PID_2, X_PID_2] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC_2, t_2-t_1, N_Simulation, Solver);
IC_3 = X_PID_2(:,end);
Parameters.gamma_2 = gamma_2*DisturbanceFactor;
[time_PID_3, X_PID_3] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC_3, t_f-t_2, N_Simulation, Solver);
time_PID = [time_PID_1, time_PID_2(2:end) + t_1, time_PID_3(2:end) + t_2];
X_PID = [X_PID_1, X_PID_2(:,2:end), X_PID_3(:,2:end)];
y_PID = X_PID(OutputIndex,:);

%% Figure Settings
ScalingFactor = 1.25;
SS = 4;
Figure_Width = 12 * SS;
Figure_Height = 3 * SS;
FontSize = ScalingFactor*3.5 * SS;
FontSize_Small = ScalingFactor*3 * SS;
FontSize_Large = ScalingFactor*6 * SS;
LineWidth = ScalingFactor*0.65 * SS;
LineWidth_Thick = ScalingFactor*1 * SS;
LineWidth_Thin = ScalingFactor*0.5 * SS;
MarkerSize = ScalingFactor*5 * SS;
MarkerSize_Small = ScalingFactor*4 * SS;
NominalColors = [55,126,184; ... % Blue
                 228,26,28; ... % Red
                 77,175,74]/255; % Green

%% Set Figure 1
Figure1_Name = 'APIDF2_Ecoli_Simulations';
Handle_Figure1 = figure();
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [-30 + Figure_Width, 10 + 3.5*Figure_Height, Figure_Width/2 * 1.25, Figure_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];

%% Axis for Simulations
Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Position = [0.07, 0.15, 0.9, 0.79];
    Handle_Axis1.Box = 'on';
    Handle_Axis1.FontSize = FontSize;
    hold(Handle_Axis1, 'on');
    grid(Handle_Axis1, 'on');
    Handle_Axis1.XMinorGrid = 'off';
    Handle_Axis1.YMinorGrid = 'off';
    Handle_Axis1.XLim = [0, t_f];
    Handle_Axis1.XTick = 0 : 4*60 : t_f;
    Handle_Axis1.XTickLabel = (0 : 4*60 : t_f)/60;
    Handle_Axis1.XLabel.String = 'Time (hr)';
    Handle_Axis1.YLabel.String = 'Output Concentration (nM)';
    Handle_Axis1.YLim = [0, 3.5*r];
plot(Handle_Axis1, time_I, y_I, 'LineWidth', LineWidth, 'Color', NominalColors(1,:));
plot(Handle_Axis1, time_PI, y_PI, 'LineWidth', LineWidth, 'Color', NominalColors(2,:));
plot(Handle_Axis1, time_PID, y_PID, 'LineWidth', LineWidth, 'Color', NominalColors(3,:));
plot(Handle_Axis1, [0, t_1, t_1, t_f], [r/TrackingFactor, r/TrackingFactor, r, r], 'LineWidth', LineWidth_Thin, 'Color', 'k', 'LineStyle', '--');
plot(Handle_Axis1, [t_1, t_1], Handle_Axis1.YLim, 'LineWidth', LineWidth_Thin, 'Color', 0.5*[1, 1, 1], 'LineStyle', '--');
plot(Handle_Axis1, [t_2, t_2], Handle_Axis1.YLim, 'LineWidth', LineWidth_Thin, 'Color', 0.5*[1, 1, 1], 'LineStyle', '--');
Handle_Legend = legend(Handle_Axis1, 'I', 'PI', 'PID', 'Set-point');
    Handle_Legend.FontSize = FontSize;
    Handle_Legend.Position = [0.1627 0.6618 0.1340 0.2632];

annotation(Handle_Figure1,'arrow',[0.37 0.37],...
    [0.964705882352941 0.941176470588236],...
    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'HeadWidth',20,...
    'HeadStyle','plain',...
    'HeadLength',13);
annotation(Handle_Figure1,'arrow',[0.670822561692126 0.670822561692126],...
    [0.964705882352941 0.941176470588236],...
    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'HeadWidth',20,...
    'HeadStyle','plain',...
    'HeadLength',13);

annotation(Handle_Figure1,'textbox',...
    [0.380728554641598 0.922055555555556 0.195945945945946 0.0823148148148148],...
    'String',{'Set-point (\mu) doubled'},...
    'FontSize',15,...
    'EdgeColor','none');
annotation(Handle_Figure1,'textbox',...
    [0.686251468860164 0.922055555555555 0.151292596944771 0.0823148148148148],...
    'String',{'Disturbance (\Delta)'},...
    'FontSize',15,...
    'EdgeColor','none');

%% Save Figures
if Save_Flag == 1  
    Handle_Figure1.Color = 'none';
    set(Handle_Figure1, 'InvertHardCopy', 'off');
  	print(Handle_Figure1, Figure1_Name, '-dpdf', '-painters');
end

