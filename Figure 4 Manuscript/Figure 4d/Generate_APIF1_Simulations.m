%% Clear Workspace
close all;
clear;
clc;

Load = 0;
if Load == 0
%% Fixed Parameters
RunTime = tic;
mu = 10;
theta = 2; 
eta = 100;
k = 3;
kappa = 1;
kappa_1 = 0.01;
n = 1;
k_1 = 2;
gamma_1 = 2;
gamma_2 = 2;
r = mu/theta;
Threshold = 0.02;

%% Swept Parameters
% Number of Steps
N = 100;  % In paper: N = 200;
% Axis Limit
alpha_max = 90;  
kappa_min = 1/0.7;
delta_max = 6;
Overshoot_max = 4;  
alpha_H = (gamma_1 * gamma_2 / k_1) * r * (1 + (r/kappa)^n);
alpha_vector = linspace(0, alpha_H, N);
kappa_vector = logspace(log10(kappa_min), 4, N);
delta_vector = linspace(0, delta_max, N);

%% Simulation Settings
IC = zeros(4,1);
N_Simulation = 1000;
Solver = 'ODE23s';
tf = 100;

%% Nominal Simulations
% Standalone AIF
clear time_vector x Parameters
Parameters.alpha = 0;
Parameters.mu = mu;
Parameters.theta = theta;
Parameters.eta = eta;
Parameters.k = k;
Parameters.n = n;
Parameters.kappa = kappa;
Parameters.beta = 0;
Parameters.gamma_c = 0;
Parameters.gamma_1 = gamma_1;
Parameters.gamma_2 = gamma_2;
Parameters.k_1 = k_1;
[time_vector, x] = DSA(StoichiometryMatrix_APIDF_Additive_2LinearSpecies(), @PropensityFunction_APIDF_Additive_2LinearSpecies, Parameters, IC, tf, N_Simulation, Solver);
x_AIF = x(2,:);
SettlingTime_AIF = Compute_SettlingTime(time_vector, x_AIF, r, Threshold);
Overshoot_AIF = max(max(x_AIF) - r, 0);
clear x;

% Additive
alpha_Nominal = 50;
clear time_vector x Parameters
Parameters.alpha = alpha_Nominal;
Parameters.mu = mu;
Parameters.theta = theta;
Parameters.eta = eta;
Parameters.k = k;
Parameters.n = n;
Parameters.kappa = kappa;
Parameters.beta = 0;
Parameters.gamma_c = 0;
Parameters.gamma_1 = gamma_1;
Parameters.gamma_2 = gamma_2;
Parameters.k_1 = k_1;
[time_vector, x] = DSA(StoichiometryMatrix_APIDF_Additive_2LinearSpecies(), @PropensityFunction_APIDF_Additive_2LinearSpecies, Parameters, IC, tf, N_Simulation, Solver);
x_Additive = x(2,:);
SettlingTime_APIF1Additive_Nominal = Compute_SettlingTime(time_vector, x_Additive, r, Threshold);
Overshoot_APIF1Additive_Nominal = max(max(x_Additive) - r, 0);
clear x;

% Multiplicative
kappa_Nominal = 1/0.43;
clear time_vector x Parameters
Parameters.kappa = kappa_Nominal;
Parameters.mu = mu;
Parameters.theta = theta;
Parameters.eta = eta;
Parameters.k = k;
Parameters.n = n;
Parameters.beta = 0;
Parameters.gamma_c = 0;
Parameters.gamma_1 = gamma_1;
Parameters.gamma_2 = gamma_2;
Parameters.k_1 = k_1;
[time_vector, x] = DSA(StoichiometryMatrix_APIDF_Multiplicative_2LinearSpecies(), @PropensityFunction_APIDF_Multiplicative_2LinearSpecies, Parameters, IC, tf, N_Simulation, Solver);
x_Multiplicative = x(2,:);
SettlingTime_APIF1Multiplicative_Nominal = Compute_SettlingTime(time_vector, x_Multiplicative, r, Threshold);
Overshoot_APIF1Multiplicative_Nominal = max(max(x_Multiplicative) - r, 0);
clear x;

% Degradation
delta_Nominal = 4.5;
clear time_vector x Parameters
Parameters.delta = delta_Nominal;
Parameters.kappa_1 = kappa_1;
Parameters.mu = mu;
Parameters.theta = theta;
Parameters.eta = eta;
Parameters.k = k;
Parameters.n = n;
Parameters.beta = 0;
Parameters.gamma_c = 0;
Parameters.gamma_1 = gamma_1;
Parameters.gamma_2 = gamma_2;
Parameters.k_1 = k_1;
[time_vector, x] = DSA(StoichiometryMatrix_APIDF_Degradation_2LinearSpecies(), @PropensityFunction_APIDF_Degradation_2LinearSpecies, Parameters, IC, tf, N_Simulation, Solver);
x_Degradation = x(2,:);
SettlingTime_APIF1Degradation_Nominal = Compute_SettlingTime(time_vector, x_Degradation, r, Threshold);
Overshoot_APIF1Degradation_Nominal = max(max(x_Degradation) - r, 0);
clear x;

%% Simulate APIF Class 1 with Additive Inhibition
clear time_vector x Parameters
Parameters.mu = mu;
Parameters.theta = theta;
Parameters.eta = eta;
Parameters.k = k;
Parameters.n = n;
Parameters.kappa = kappa;
Parameters.beta = 0;
Parameters.gamma_c = 0;
Parameters.gamma_1 = gamma_1;
Parameters.gamma_2 = gamma_2;
Parameters.k_1 = k_1;
SettlingTime_APIF1Additive = zeros(N,1);
Overshoot_APIF1Additive = zeros(N,1);
for i = 1 : N
    i
    Parameters.alpha = alpha_vector(i);
    [time_vector, x] = DSA(StoichiometryMatrix_APIDF_Additive_2LinearSpecies(), @PropensityFunction_APIDF_Additive_2LinearSpecies, Parameters, IC, tf, N_Simulation, Solver);
    SettlingTime_APIF1Additive(i) = Compute_SettlingTime(time_vector, x(2,:), r, Threshold);
    Overshoot_APIF1Additive(i) = max(max(x(2,:)) - r, 0);
end

%% Simulate APIF Class 1 with Multiplicative Inhibition
clear time_vector x Parameters
Parameters.mu = mu;
Parameters.theta = theta;
Parameters.eta = eta;
Parameters.k = k;
Parameters.n = n;
Parameters.beta = 0;
Parameters.gamma_c = 0;
Parameters.gamma_1 = gamma_1;
Parameters.gamma_2 = gamma_2;
Parameters.k_1 = k_1;
SettlingTime_APIF1Multiplicative = zeros(N,1);
Overshoot_APIF1Multiplicative = zeros(N,1);
for i = 1 : N
    i
    Parameters.kappa = kappa_vector(i);
    [time_vector, x] = DSA(StoichiometryMatrix_APIDF_Multiplicative_2LinearSpecies(), @PropensityFunction_APIDF_Multiplicative_2LinearSpecies, Parameters, IC, tf, N_Simulation, Solver);
    SettlingTime_APIF1Multiplicative(i) = Compute_SettlingTime(time_vector, x(2,:), r, Threshold);
    Overshoot_APIF1Multiplicative(i) = max(max(x(2,:)) - r, 0);
end

%% Simulate APIF Class 1 with Degradation Inhibition
clear time_vector x Parameters
Parameters.kappa_1 = kappa_1;
Parameters.mu = mu;
Parameters.theta = theta;
Parameters.eta = eta;
Parameters.k = k;
Parameters.n = n;
Parameters.beta = 0;
Parameters.gamma_c = 0;
Parameters.gamma_1 = gamma_1;
Parameters.gamma_2 = gamma_2;
Parameters.k_1 = k_1;
SettlingTime_APIF1Degradation = zeros(N,1);
Overshoot_APIF1Degradation = zeros(N,1);
for i = 1 : N
    i
    Parameters.delta = delta_vector(i);
    [time_vector, x] = DSA(StoichiometryMatrix_APIDF_Degradation_2LinearSpecies(), @PropensityFunction_APIDF_Degradation_2LinearSpecies, Parameters, IC, tf, N_Simulation, Solver);
    SettlingTime_APIF1Degradation(i) = Compute_SettlingTime(time_vector, x(2,:), r, Threshold);
    Overshoot_APIF1Degradation(i) = max(max(x(2,:)) - r, 0);
end

RunTime = toc(RunTime);
save APIF1_DeterministicSimulations

else
load APIF1_DeterministicSimulations
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
Figure1_Name = 'APIF1_Settling_Overshoot';
Handle_Figure1 = figure();
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [-30, 50, Figure_Width, Figure_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];
    
%% Set Figure 2
Figure2_Name = 'APIF1_NominalSimulations';
Handle_Figure2 = figure();
    Handle_Figure2.PaperUnits = 'centimeters';
    Handle_Figure2.Units = 'centimeters';
    Handle_Figure2.Position = [-30 + Figure_Width, 50, Figure_Width, Figure_Height];
    Handle_Figure2.PaperPositionMode = 'auto';
    Handle_Figure2.PaperSize = [Handle_Figure2.PaperPosition(3), Handle_Figure2.PaperPosition(4)];
    
%% Axis for APIF1 (Additive) Settling Time
Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Position = [0.09, 0.18, 0.84, 0.73];
    Handle_Axis1.Box = 'off';
    Handle_Axis1.FontSize = FontSize;
    hold(Handle_Axis1, 'on');
    grid(Handle_Axis1, 'off');
    Handle_Axis1.YGrid = 'on';
    Handle_Axis1.XMinorGrid = 'off';
    Handle_Axis1.YMinorGrid = 'off';
    Handle_Axis1.XLabel.String = '$\alpha$';
    Handle_Axis1.XLabel.Interpreter = 'latex';
    Handle_Axis1.YLabel.String = 'Settling Time';
    Handle_Axis1.XScale = 'linear';
    Handle_Axis1.YScale = 'linear'; 
    Handle_Axis1.XColor = MyBlue;
    Handle_Axis1.XLim = [0, alpha_max];
    Handle_Axis1.XTick = [0, 30, 60, 90];
    Handle_Axis1.XLabel.Position(1:2) = [alpha_max/2, 0];
plot(Handle_Axis1, alpha_vector, SettlingTime_APIF1Additive, 'LineWidth', LineWidth, 'Color', MyBlue);
plot(Handle_Axis1, alpha_H*[1,1], Handle_Axis1.YLim, 'LineWidth', LineWidth/2, 'LineStyle', ':', 'Color', MyBlue);
plot(Handle_Axis1, alpha_Nominal, SettlingTime_APIF1Additive_Nominal, 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', MyBlue, 'MarkerFaceColor', MyBlue);
Handle_Annotation1 = annotation(Handle_Figure1, 'textbox', 'LineStyle', 'none');
    Handle_Annotation1.Interpreter = 'latex';
    Handle_Annotation1.String = 'Unstable';
    Handle_Annotation1.FontName = 'Helvetica';
    Handle_Annotation1.FontSize = FontSize/1.2;
    Handle_Annotation1.Position = [0.66 0.55 0 0];
    Handle_Annotation1.Color = MyBlue;
Handle_Annotation2 = annotation(Handle_Figure1, 'arrow', 'LineStyle', '-', 'Color', MyBlue);
    Handle_Annotation2.Position = [0.66, 0.49, 0.18, 0];

%% Axis for APIF1 (Multiplicative) Settling Time
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
    Handle_Axis2.XTick = [0, 0.25, 0.5];
    Handle_Axis2.XLabel.Position(1:2) = [(1/kappa_min)/2, Handle_Axis1.YLim(2)];
plot(Handle_Axis2, 1./kappa_vector, SettlingTime_APIF1Multiplicative, 'LineWidth', LineWidth, 'Color', MyRed);
plot(Handle_Axis2, 1/kappa_Nominal, SettlingTime_APIF1Multiplicative_Nominal, 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', MyRed, 'MarkerFaceColor', MyRed);

%% Axis for APIF1 (Degradation) Settling Time
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
plot(Handle_Axis3, delta_vector, SettlingTime_APIF1Degradation, 'LineWidth', LineWidth, 'Color', MyGreen);
plot(Handle_Axis3, delta_Nominal, SettlingTime_APIF1Degradation_Nominal, 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', MyGreen, 'MarkerFaceColor', MyGreen);

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
    Handle_Axis3_Prime.XLabel.String = '$\delta$';
    Handle_Axis3_Prime.XLabel.Interpreter = 'latex';
    Handle_Axis3_Prime.XLim = Handle_Axis3.XLim;
    Handle_Axis3_Prime.XTick = linspace(delta_vector(1), delta_vector(end), 4);
    Handle_Axis3_Prime.XLabel.Position(1:2) = [delta_max/2, 0];
    
%% Axis for APIF1 (Additive) Overshoot
Handle_Axis4 = axes(Handle_Figure1);
    Handle_Axis4.Position = Handle_Axis1.Position;
    Handle_Axis4.Box = 'off';
    Handle_Axis4.FontSize = FontSize;
    hold(Handle_Axis4, 'on');
    grid(Handle_Axis4, 'off');
    Handle_Axis4.XMinorGrid = 'off';
    Handle_Axis4.YMinorGrid = 'off';
    Handle_Axis4.XScale = 'linear';
    Handle_Axis4.YScale = 'linear'; 
    Handle_Axis4.Color = 'none';
    Handle_Axis4.XColor = 'none';
    Handle_Axis4.YAxisLocation = 'right';
    Handle_Axis4.XLim = [0, alpha_max];
    Handle_Axis4.YLim = [0, Overshoot_max];
    Handle_Axis4.YLabel.String = 'Overshoot';
plot(Handle_Axis4, alpha_vector, Overshoot_APIF1Additive, 'LineWidth', LineWidth, 'Color', MyBlue, 'LineStyle', '--');
plot(Handle_Axis4, alpha_Nominal, Overshoot_APIF1Additive_Nominal, 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', MyBlue, 'MarkerFaceColor', MyBlue);

%% Axis for APIF1 (Multiplicative) Overshoot
Handle_Axis5 = axes(Handle_Figure1);
    Handle_Axis5.Position = Handle_Axis1.Position;
    Handle_Axis5.Box = 'off';
    Handle_Axis5.FontSize = FontSize;
    hold(Handle_Axis5, 'on');
    grid(Handle_Axis5, 'off');
    Handle_Axis5.XMinorGrid = 'off';
    Handle_Axis5.YMinorGrid = 'off';
    Handle_Axis5.XScale = 'linear';
    Handle_Axis5.YScale = 'linear'; 
    Handle_Axis5.Color = 'none';
    Handle_Axis5.XColor = 'none';
    Handle_Axis5.YAxisLocation = 'right';
    Handle_Axis5.XLim = [0, 1/kappa_min];
    Handle_Axis5.YLim = [0, Overshoot_max];
plot(Handle_Axis5, 1./kappa_vector, Overshoot_APIF1Multiplicative, 'LineWidth', LineWidth, 'Color', MyRed, 'LineStyle', '--');
plot(Handle_Axis5, 1./kappa_Nominal, Overshoot_APIF1Multiplicative_Nominal, 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', MyRed, 'MarkerFaceColor', MyRed);

%% Axis for APIF1 (Degradation) Overshoot
Handle_Axis6 = axes(Handle_Figure1);
    Handle_Axis6.Position = Handle_Axis1.Position;
    Handle_Axis6.Box = 'off';
    Handle_Axis6.FontSize = FontSize;
    hold(Handle_Axis6, 'on');
    grid(Handle_Axis6, 'off');
    Handle_Axis6.XMinorGrid = 'off';
    Handle_Axis6.YMinorGrid = 'off';
    Handle_Axis6.XScale = 'linear';
    Handle_Axis6.YScale = 'linear'; 
    Handle_Axis6.Color = 'none';
    Handle_Axis6.XColor = 'none';
    Handle_Axis6.YAxisLocation = 'right';
    Handle_Axis6.XLim = [0, delta_max];
    Handle_Axis6.YLim = [0, Overshoot_max];
plot(Handle_Axis6, delta_vector, Overshoot_APIF1Degradation, 'LineWidth', LineWidth, 'Color', MyGreen, 'LineStyle', '--');
plot(Handle_Axis6, delta_Nominal, Overshoot_APIF1Degradation_Nominal, 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', MyGreen, 'MarkerFaceColor', MyGreen);

%% Legend
Handle_Axis_Legend1 = axes(Handle_Figure1);
    Handle_Axis_Legend1.Position = Handle_Axis1.Position;
    hold(Handle_Axis_Legend1, 'on');
patch(Handle_Axis_Legend1, -1, -1, MyBlue, 'EdgeColor' , 'none');
patch(Handle_Axis_Legend1, -1, -1, MyRed, 'EdgeColor' , 'none');
patch(Handle_Axis_Legend1, -1, -1, MyGreen, 'EdgeColor' , 'none');
Handle_Legend1 = legend(Handle_Axis_Legend1, {'Additive PI', 'Multiplicative PI', 'Degradation PI'});
    Handle_Legend1.FontSize = FontSize/1.2;
    Handle_Legend1.Position = [0.34, 0.82, 0, 0];
    Handle_Axis_Legend1.Visible = 'off';
    
Handle_Axis_Legend2 = axes(Handle_Figure1);
    Handle_Axis_Legend2.Position = Handle_Axis1.Position;
    hold(Handle_Axis_Legend2, 'on');
plot(Handle_Axis_Legend2, -1, -1, 'LineWidth', LineWidth, 'Color', [0.5, 0.5, 0.5]);
plot(Handle_Axis_Legend2, -1, -1, 'LineWidth', LineWidth, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--');
Handle_Legend2 = legend(Handle_Axis_Legend2, {'Settling Time', 'Overshoot'});
    Handle_Legend2.FontSize = FontSize/1.2;
    Handle_Legend2.Position = [0.6, 0.82, 0, 0];
    Handle_Axis_Legend2.Visible = 'off';
    
%% Axis for Nominal Simulations
Handle_Axis7 = axes(Handle_Figure2);
    Handle_Axis7.Position = [0.09, 0.18, 0.85, 0.73];
    Handle_Axis7.Box = 'on';
    Handle_Axis7.FontSize = FontSize;
    hold(Handle_Axis7, 'on');
    grid(Handle_Axis7, 'on');
    Handle_Axis7.XMinorGrid = 'off';
    Handle_Axis7.YMinorGrid = 'off';
    Handle_Axis7.XLabel.String = '$t$';
    Handle_Axis7.XLabel.Interpreter = 'latex';
    Handle_Axis7.YLabel.String = '$x_2(t)$';
    Handle_Axis7.YLabel.Interpreter = 'latex';
    Handle_Axis7.XScale = 'linear';
    Handle_Axis7.YScale = 'linear'; 
    Handle_Axis7.XLim = [0, 50];
    Handle_Axis7.XLabel.Position(1:2) = [tf/4, 0];
plot(Handle_Axis7, time_vector, x_AIF, 'LineWidth', LineWidth/2, 'Color', 'k');
plot(Handle_Axis7, time_vector, x_Additive, 'LineWidth', LineWidth, 'Color', MyBlue);
plot(Handle_Axis7, time_vector, x_Multiplicative, 'LineWidth', LineWidth, 'Color', MyRed);
plot(Handle_Axis7, time_vector, x_Degradation, 'LineWidth', LineWidth, 'Color', MyGreen, 'LineStyle', '--');
Handle_Legend3 = legend(Handle_Axis7, {'Standalone I', 'Additive PI', 'Multiplicative PI', 'Degradation PI'});
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

end
