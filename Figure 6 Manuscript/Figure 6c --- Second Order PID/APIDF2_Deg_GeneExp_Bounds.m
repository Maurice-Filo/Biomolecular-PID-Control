%% Clear Workspace
close all
clear
clc
Save_Flag = 0;

%% Functions
StoichiometryMatrix = StoichiometryMatrix_APIDF2_Deg_GeneExp();
PropensityFunction = @PropensityFunction_APIDF2_Deg_GeneExp;
FixedPoint = @APIDF2_Deg_FixedPoint;
Params2Gains = @APIDF2_Deg_Params2Gains;
Gains2Params = @APIDF2_Deg_Gains2Params;
SupportingInput = @SupportingInput_GeneExp;
TransferFunction = @APIDF2_Deg_TF;
OutputIndex = 2;

%% Plant Parameters
Parameters.gamma_1 = 0.1;
Parameters.gamma_2 = 0.1;
Parameters.k_1 = 0.1;

%% Controller Parameters
Parameters.mu = 10;
Parameters.kappa = 1e-5;
Parameters.n = 1;
Parameters.delta_c = 0;
r = 5;
K_S = Parameters.mu / r;

%% Disturbance
Tracking_Factor = 2;
r_New = r * Tracking_Factor;
Parameters_New = Parameters;
Parameters_New.mu = Parameters.mu * Tracking_Factor;
u_bar = SupportingInput(Parameters, r_New);

%% Simulation Settings
N_Simulation = 1000;
Solver = 'ODE15s';
t_1 = 100;
t_f = 130;

%% Sweep the Pole
N = 1000;
N_Nominal = 4;
gamma_bar = (Parameters.gamma_1 + Parameters.gamma_2) / 2;
a_Sweep = linspace(gamma_bar, gamma_bar * (sqrt(2) + 2), N);
Indeces_Nominal = round(linspace(1, N, N_Nominal));

%% Compute Gains and Biological Parameters for each place pole
Gains = cell(N,1);
omega_c_vector = zeros(N,1);
K_P_vector = zeros(N,1);
K_I_vector = zeros(N,1);
K_D_vector = zeros(N,1);
CoverageFlag = zeros(N,1);
BioParameters = cell(N,1);
Feasibility_Flag = zeros(N,1);
eta_vector = zeros(N,1);
delta_vector = zeros(N,1);
k_vector = zeros(N,1);
beta_vector = zeros(N,1);
for i = 1 : N
    a = a_Sweep(i);
    omega_c = 4*a - (Parameters.gamma_1 + Parameters.gamma_2);
    K_I = a^4 / (Parameters.k_1 * K_S * omega_c);
    K_P = 4*a^3 / Parameters.k_1 / omega_c - Parameters.gamma_1*Parameters.gamma_2/Parameters.k_1;
    K_D = (6*a^2 - Parameters.gamma_1*Parameters.gamma_2) / Parameters.k_1 / omega_c - (Parameters.gamma_1 + Parameters.gamma_2) / Parameters.k_1;
    Gains{i}.K_P = K_P;
    Gains{i}.K_I = K_I;
    Gains{i}.K_D = K_D; 
    Gains{i}.omega_c = omega_c;
    Gains{i}.K_F = 0.1;
    Gains{i}.K_S = K_S; 
    omega_c_vector(i) = omega_c;
    K_I_vector(i) = K_I;
    K_P_vector(i) = K_P;
    K_D_vector(i) = K_D;
    CoverageFlag(i) = (K_P < K_D*omega_c) && (K_I < (omega_c/Parameters_New.mu/4) * (u_bar + r_New*K_P)^2 / (u_bar + r_New*K_D*omega_c));
    [BioParameters{i}, Feasibility_Flag(i)] = Gains2Params(Gains{i}, Parameters_New, r_New, SupportingInput);
    eta_vector(i) = BioParameters{i}.eta;
    delta_vector(i) = BioParameters{i}.delta;
    k_vector(i) = BioParameters{i}.k;
    beta_vector(i) = BioParameters{i}.beta;
end

%% Simulations
a_Nominal = zeros(N_Nominal,1);
Gains_Nominal = cell(N_Nominal,1);
BioParameters_Nominal = cell(N_Nominal,1);
x_APIDF2 = zeros(N_Nominal, 2*N_Simulation-1);
for i = 1 : N_Nominal
    Index = Indeces_Nominal(i);
    a_Nominal(i) = a_Sweep(Index);
    Parameters = BioParameters{Index};
    Parameters.mu = Parameters.mu / Tracking_Factor;
    IC = FixedPoint(Parameters, SupportingInput);
    [time_vector_1, X_1] = DSA(StoichiometryMatrix, PropensityFunction, Parameters, IC, t_1, N_Simulation, Solver);
    [time_vector_2, X_2] = DSA(StoichiometryMatrix, PropensityFunction, BioParameters{Index}, X_1(:,end), t_f-t_1, N_Simulation, Solver);
    X = [X_1, X_2(:,2:end)];
    time_vector = [time_vector_1, t_1 + time_vector_2(2:end)];
    x_APIDF2(i,:) = X(OutputIndex,:);
    Gains_Nominal{i} = Params2Gains(BioParameters{Index}, SupportingInput);
    BioParameters_Nominal{i} = BioParameters{Index};
end

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
LineWidth_Thin = ScalingFactor*0.5 * SS;
MarkerSize = ScalingFactor*5 * SS;
MarkerSize_Small = ScalingFactor*4 * SS;
NominalColors = [55,126,184; ... % Blue
                 77,175,74; ... % Green
                 255,127,0; ... % Orange
                 228,26,28]/255; ... % Red
             
             
%% Set Figure 1
Figure1_Name = 'APIDF2_Deg_GeneExp_Poles';
Handle_Figure1 = figure();
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [-40, 10 + 3.5*Figure_Height, Figure_Width/2.5, Figure_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];
    
%% Set Figure 2
Figure2_Name = 'APIDF2_Deg_GeneExp_Gains';
Handle_Figure2 = figure();
    Handle_Figure2.PaperUnits = 'centimeters';
    Handle_Figure2.Units = 'centimeters';
    Handle_Figure2.Position = [-40 + Figure_Width/2.5, 10 + 3.5*Figure_Height, Figure_Width/2.5, Figure_Height];
    Handle_Figure2.PaperPositionMode = 'auto';
    Handle_Figure2.PaperSize = [Handle_Figure2.PaperPosition(3), Handle_Figure2.PaperPosition(4)];
    
%% Set Figure 3
Figure3_Name = 'APIDF2_Deg_GeneExp_BioParameters';
Handle_Figure3 = figure();
    Handle_Figure3.PaperUnits = 'centimeters';
    Handle_Figure3.Units = 'centimeters';
    Handle_Figure3.Position = [-40 + 2*Figure_Width/2.5, 10 + 3.5*Figure_Height, Figure_Width/3.2 * 1.25, Figure_Height];
    Handle_Figure3.PaperPositionMode = 'auto';
    Handle_Figure3.PaperSize = [Handle_Figure3.PaperPosition(3), Handle_Figure3.PaperPosition(4)];

%% Set Figure 4
Figure4_Name = 'APIDF2_Deg_GeneExp_Simulations';
Handle_Figure4 = figure();
    Handle_Figure4.PaperUnits = 'centimeters';
    Handle_Figure4.Units = 'centimeters';
    Handle_Figure4.Position = [-30 + Figure_Width, 10 + 3.5*Figure_Height, Figure_Width/2 * 1.25, Figure_Height];
    Handle_Figure4.PaperPositionMode = 'auto';
    Handle_Figure4.PaperSize = [Handle_Figure4.PaperPosition(3), Handle_Figure4.PaperPosition(4)];

%% Axis for Poles
Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Position = [0.02, 0.16, 0.78, 0.76];
    Handle_Axis1.Box = 'on';
    Handle_Axis1.BoxStyle = 'full';
    Handle_Axis1.LineWidth = LineWidth_Thin;
    Handle_Axis1.FontSize = FontSize_Small*1.25;
    hold(Handle_Axis1, 'on');
    grid(Handle_Axis1, 'on');
    Handle_Axis1.XMinorGrid = 'off';
    Handle_Axis1.YMinorGrid = 'off';
    Handle_Axis1.XLim = [-0.8,0.1];
    Handle_Axis1.XTick = [-0.6, -round(gamma_bar*(2 + sqrt(2)),2), -gamma_bar, 0];
    Handle_Axis1.YLim = [-0.4, 0.4];
    Handle_Axis1.YTick = [-0.2, 0, 0.2];
    Handle_Axis1.XLabel.String = 'Real Axis';
    Handle_Axis1.XLabel.FontSize = FontSize;
    Handle_Axis1.YLabel.String = 'Imaginary Axis';
    Handle_Axis1.YLabel.FontSize = FontSize;
    colormap(Handle_Axis1, 'turbo');
    Handle_Axis1.YAxisLocation = 'Right';
    Handle_Axis1.Title.String = 'Pole Placement';
    Handle_Axis1.Title.FontSize = FontSize;
scatter(Handle_Axis1, -a_Sweep, 0*a_Sweep, MarkerSize*3, a_Sweep, 'filled');
for i = 1 : N_Nominal
    plot(Handle_Axis1, -a_Nominal(i), 0, 'Marker', 'o', 'MarkerSize', MarkerSize_Small, 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', NominalColors(i,:));
end
plot(Handle_Axis1, [0, 0], Handle_Axis1.YLim, 'LineWidth', LineWidth_Thin, 'Color', 'k', 'LineStyle', ':');
plot(Handle_Axis1, -gamma_bar * [1, 1], [Handle_Axis1.YLim(1), 0], 'LineWidth', LineWidth_Thin, 'Color', 'k', 'LineStyle', '--');
plot(Handle_Axis1, -gamma_bar * (2 + sqrt(2)) * [1, 1], [Handle_Axis1.YLim(1), 0], 'LineWidth', LineWidth_Thin, 'Color', 'k', 'LineStyle', '--');
Handle_Annotation1 = annotation(Handle_Figure1, 'textarrow', [0.258680555555556,0.413020833333334], [0.640833333333333,0.163652392947103]);
    Handle_Annotation1.String = '$-(2 + \sqrt{2}) \frac{\gamma_1 + \gamma_2}{2}$';
    Handle_Annotation1.Interpreter = 'latex';
    Handle_Annotation1.FontSize = FontSize;
    Handle_Annotation1.HorizontalAlignment = 'center';
    Handle_Annotation1.Color = 'r';
    Handle_Annotation1.TextColor = 'k';
Handle_Annotation2 = annotation(Handle_Figure1, 'textarrow', [0.572873263888889,0.625238715277779], [0.787824074074074,0.163652392947103]);
    Handle_Annotation2.String = '$-\frac{\gamma_1 + \gamma_2}{2}$';
    Handle_Annotation2.Interpreter = 'latex';
    Handle_Annotation2.FontSize = FontSize;
    Handle_Annotation2.HorizontalAlignment = 'center';
    Handle_Annotation2.Color = 'b';
    Handle_Annotation2.TextColor = 'k';
    
%% Axis for Gains
Handle_Axis2 = axes(Handle_Figure2);
    Handle_Axis2.Position = [0.17, 0.1, 0.7, 0.82];
    Handle_Axis2.Box = 'on';
    Handle_Axis2.BoxStyle = 'full';
    Handle_Axis2.LineWidth = LineWidth_Thin;
    Handle_Axis2.FontSize = FontSize;
    hold(Handle_Axis2, 'on');
    grid(Handle_Axis2, 'on');
    Handle_Axis2.XMinorGrid = 'off';
    Handle_Axis2.YMinorGrid = 'off';
    Handle_Axis2.XLabel.String = '$K_P$';
    Handle_Axis2.YLabel.String = '$K_I$';
    Handle_Axis2.ZLabel.String = '$K_D$';
    Handle_Axis2.XLabel.Interpreter = 'latex';
    Handle_Axis2.YLabel.Interpreter = 'latex';
    Handle_Axis2.ZLabel.Interpreter = 'latex';
    view(Handle_Axis2, [40, 40])
    colormap(Handle_Axis2, 'turbo');
 	Handle_Axis2.Title.String = 'PID Gains';
scatter3(Handle_Axis2, K_P_vector, K_I_vector, K_D_vector, MarkerSize*3, linspace(0, 1, N), 'filled');
for i = 1 : N_Nominal
    plot3(Handle_Axis2, Gains_Nominal{i}.K_P, Gains_Nominal{i}.K_I, Gains_Nominal{i}.K_D, 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', NominalColors(i,:));
end

%% Axis for BioParameters
Handle_Axis3 = axes(Handle_Figure3);
    Handle_Axis3.Position = [0.17, 0.1, 0.7, 0.82];
    Handle_Axis3.Box = 'on';
    Handle_Axis3.BoxStyle = 'full';
    Handle_Axis3.LineWidth = LineWidth_Thin;
    Handle_Axis3.FontSize = FontSize;
    hold(Handle_Axis3, 'on');
    grid(Handle_Axis3, 'on');
    Handle_Axis3.XMinorGrid = 'off';
    Handle_Axis3.YMinorGrid = 'off';
    Handle_Axis3.XLabel.String = '$\delta$';
    Handle_Axis3.YLabel.String = '$k$';
    Handle_Axis3.ZLabel.String = '$\beta$';
    Handle_Axis3.XLabel.Interpreter = 'latex';
    Handle_Axis3.YLabel.Interpreter = 'latex';
    Handle_Axis3.ZLabel.Interpreter = 'latex';
    view(Handle_Axis3, [40, 40])
    colormap(Handle_Axis3, 'turbo');
    Handle_Axis3.XTick = [0, 2, 4];
    Handle_Axis3.Title.String = 'Biomolecular Parameters';
scatter3(Handle_Axis3, delta_vector, k_vector, beta_vector, MarkerSize*3, linspace(0, 1, N), 'filled');
for i = 1 : N_Nominal
    plot3(Handle_Axis3, BioParameters_Nominal{i}.delta, BioParameters_Nominal{i}.k, BioParameters_Nominal{i}.beta, 'Marker', 'o', 'MarkerSize', MarkerSize, 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', NominalColors(i,:));
end

%% Axis for Simulations
Handle_Axis4 = axes(Handle_Figure4);
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
    Handle_Axis4.YLim = [0.9*r, 1.1*r_New];
    Handle_Axis4.Title.String = 'Simulations';
plot(Handle_Axis4, [0, t_1, t_1, t_f], [r, r, r_New, r_New], 'LineWidth', LineWidth_Thin, 'Color', 'k', 'LineStyle', '--');
for i = 1 : N_Nominal
    plot(Handle_Axis4, time_vector, x_APIDF2(i,:), 'LineWidth', LineWidth, 'Color', NominalColors(i,:));
end

%% Save Figures
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

