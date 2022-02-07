%% Clear Workspace
close all;
clear;
clc;
Colors = lines(10);
MyBlue = [0, 130, 200]/255; 
MyRed = [230, 25, 75]/255; 
MyGreen = [60, 180, 75]/255;
MyCyan = 0.4*[70, 240, 240]/255 + 0.6*[1, 1, 1];
MyLime = 0.4*[210. 245, 60]/255 + 0.6*[1, 1, 1];
SS = 4; % Screen Scale
Save_Flag = 0;

addpath('Functions');

%% Plant Parameters
PlantParameters.Ingestion.k1 = 0;
PlantParameters.ICSmax = 100;
PlantParameters.ICt = 3120;
PlantParameters.IS = 5;
PlantParameters.BileSaltRelease.k1 = 6;
PlantParameters.BileSaltReturn.k1 = 4.29;
PlantParameters.BileSaltExcretion.k1 = 0.856;
PlantParameters.k5 = 2.66;
PlantParameters.k6 = 0.6545*0.0005286;
PlantParameters.k7 = 0.0005286;
PlantParameters.k8 = 0.0005;
PlantParameters.BCRmax = 2000;
PlantParameters.BCRt = 55326;
PlantParameters.BS = 5;
PlantParameters.HCSmax = 500;
PlantParameters.HCSt = 93925;
PlantParameters.HS = 5;
PlantParameters.k9 = 1;
PlantParameters.k10 = 5.998;
PlantParameters.k11 = 0.005;
PlantParameters.VLDLCholesterolFormation.k1 = 0.016;
PlantParameters.khrs = 100;
PlantParameters.HepaticLDLReceptorDegradation.k1 = 0.01;
PlantParameters.VLDLCholesterolReUptake.k1 = 0.0496;
PlantParameters.k15 = 0.43;
PlantParameters.IDLCholesterolReUptake.k1 = 0.054;
PlantParameters.k17 = 0.38;
PlantParameters.k18 = 0.068;
PlantParameters.ReceptorIndependentHepaticUptake.k1 = 0.005;
PlantParameters.k20 = 0.00675;
PlantParameters.ReceptorIndependentPeripheralUptake.k1 = 5e-06;
PlantParameters.kprs = 100;
PlantParameters.PeripheralLDLReceptorDegradation.k1 = 0.01;
PlantParameters.k23 = 0.017386;
PlantParameters.k24 = 0.1068;
PlantParameters.PeripheralSteroidProduction.k1 = 0.0005;
PlantParameters.k26 = 1.5e-05;
PlantParameters.PCSmax = 500;
PlantParameters.PPCt = 80342;
PlantParameters.PCSS = 5;
PlantParameters.k27 = 0.01;
PlantParameters.k28 = 0.001;
PlantParameters.k29 = 0.05;
PlantParameters.Intake = 1;
PlantParameters.Intestine = 1;
PlantParameters.HepaticTissue = 1;
PlantParameters.PeripheralTissue = 1;
PlantParameters.Plasma = 1;
PlantParameters.Excreted = 1;

%% Integral Controller Parameters
AIFParameters.mu = 100;
AIFParameters.theta = 1;
AIFParameters.eta = 0.1;
AIFParameters.k = 100;
AIFParameters.delta_c = 0;

%% Plant Initial Condition
IC = ...
[...
%   Initial Condition   Name        Index
    304; ...         	DC          1
    3150; ...          	IC          2
    0; ...             	ICS         3
    467; ...           	IBS         4
    400; ...            HBS         5
  	60000; ...          HFC         6
    0; ...              HCS         7
    10000; ...          HCE         8
    100; ...            ACAT        9
    100; ...            CEH         10
  	0; ...              HNHDLS      11
    100; ...            HLDLRs      12
    600; ...            HLDLRsS     13
    0; ...              HLDLRD      14
    100; ...            SRB1        15
  	57516; ...          PFC         16
    100; ...            PLDLRs      17
    575.16; ...         PLDLRsS     18
    0; ...              PLDLRD      19
    9363; ...           PCE         20
 	0; ...              PSS         21
    0; ...              PCS         22
    0; ...              INHDLS      23
    100; ...            NHDL        24
    20; ...             VLDLC       25 
	20; ...             IDLC        26
    100; ...            LPL         27
    100; ...         	LDLC        28
    100; ...            HSL         29
    45; ...             HDLC        30 
   	100; ...            LCAT        31
    100; ...            CETP        32
    0; ...              EBS         33
    0 ...               EC          34
];
Input_Index = 2;
Output_Index = 28;

%% Simulation Settings
N_Simulation = 1000;
Solver = 'ODE15s';
tf = 1600;
t0 = 900;

%% Disturbances
t1 = 1000; DisturbedParam1 = {'Ingestion'; 'k1'}; Factor1 = 1;
t2 = 1250; DisturbedParam2 = {'k6'}; Factor2 = 0.25*0.0005286;

PlantParameters1 = PlantParameters;
    if length(DisturbedParam1) == 1
        PlantParameters1.(DisturbedParam1{1}) = PlantParameters.(DisturbedParam1{1}) + Factor1;
    else
        PlantParameters1.(DisturbedParam1{1}).(DisturbedParam1{2}) = PlantParameters.(DisturbedParam1{1}).(DisturbedParam1{2}) + Factor1;
    end
PlantParameters2 = PlantParameters1;
    if length(DisturbedParam2) == 1
        PlantParameters2.(DisturbedParam2{1}) = PlantParameters1.(DisturbedParam2{1}) + Factor2;
    else
        PlantParameters2.(DisturbedParam2{1}).(DisturbedParam2{2}) = PlantParameters1.(DisturbedParam2{1}).(DisturbedParam2{2}) + Factor2;
    end
N0 = round(t1 * N_Simulation /tf);
N1 = round((t2-t1) * N_Simulation /tf);
N2 = N_Simulation - N0 - N1;

%% Simulation (Open Loop)
% Network Parameters
NetworkParameters = PlantParameters;
NetworkParameters1 = PlantParameters1;
NetworkParameters2 = PlantParameters2;
% Initial Conditions
InitialCondition = IC;
% Stoichiometry Matrix and Propensity Function
StoichiometryMatrix = StoichiometryMatrix_OL();
PropensityFunction = @PropensityFunction_OL;
% Deterministic Simulation
TStart = tic;
% Before the Disturbances
[time_vector1, x1] = DSA(StoichiometryMatrix, PropensityFunction, NetworkParameters, InitialCondition, t1, N0, Solver);
% First Disturbance at time t1
[time_vector2, x2] = DSA(StoichiometryMatrix, PropensityFunction, NetworkParameters1, x1(:,end), t2-t1, N1, Solver);
% Second Disturbance at time t2
[time_vector3, x3] = DSA(StoichiometryMatrix, PropensityFunction, NetworkParameters2, x2(:,end), tf - t2, N2, Solver);
SimTime_OL = toc(TStart);
% Construct Full Solution
time_vector_OL = [time_vector1, t1 + time_vector2(2:end), t2 + time_vector3(2:end)];
x_OL = [x1, x2(:,2:end), x3(:,2:end)];
clear x1 x2 x3 time_vector1 time_vector2 time_vector3 ...
      NetworkParameters1 NetworkParameters2 NetworkParameters ...
      StoichiometryMatrix PropensityFunction;

%% Simulation (AIF)
% Controller Parameters
ControllerParameters = AIFParameters;
% Network Parameters
NetworkParameters = CombineStructures(PlantParameters, ControllerParameters);
NetworkParameters1 = CombineStructures(PlantParameters1, ControllerParameters);
NetworkParameters2 = CombineStructures(PlantParameters2, ControllerParameters);
% Initial Condition
InitialCondition = [IC; 0; 0];
% Stoichiometry Matrix and Propensity Function
StoichiometryMatrix = StoichiometryMatrix_AIF_PositiveZ1();
PropensityFunction = @PropensityFunction_AIF_PositiveZ1;
% Deterministic Simulation
TStart = tic;
% Before the Disturbances
[time_vector1, x1] = DSA(StoichiometryMatrix, PropensityFunction, NetworkParameters, InitialCondition, t1, N0, Solver);
% First Disturbance at time t1
[time_vector2, x2] = DSA(StoichiometryMatrix, PropensityFunction, NetworkParameters1, x1(:,end), t2-t1, N1, Solver);
% Second Disturbance at time t2
[time_vector3, x3] = DSA(StoichiometryMatrix, PropensityFunction, NetworkParameters2, x2(:,end), tf - t2, N2, Solver);
SimTime_AIF = toc(TStart);
% Construct Full Solution
time_vector_AIF = [time_vector1, t1 + time_vector2(2:end), t2 + time_vector3(2:end)];
x_AIF = [x1, x2(:,2:end), x3(:,2:end)];
clear x1 x2 x3 time_vector1 time_vector2 time_vector3 ...
      NetworkParameters1 NetworkParameters2 NetworkParameters ...
      StoichiometryMatrix PropensityFunction;

%% Simulation (APIF, Degradation)
% Controller Parameters
ControllerParameters = AIFParameters;
ControllerParameters.beta = 0;
ControllerParameters.theta = ControllerParameters.theta + ControllerParameters.beta;
ControllerParameters.delta = 0.01;
ControllerParameters.n = 1;
% Network Parameters
NetworkParameters = CombineStructures(PlantParameters, ControllerParameters);
NetworkParameters1 = CombineStructures(PlantParameters1, ControllerParameters);
NetworkParameters2 = CombineStructures(PlantParameters2, ControllerParameters);
% Initial Condition
InitialCondition = [IC; 0; 0];
% Stoichiometry Matrix and Propensity Function
StoichiometryMatrix = StoichiometryMatrix_APIDF_PositiveZ1_Deg();
PropensityFunction = @PropensityFunction_APIDF_PositiveZ1_Deg;
% Deterministic Simulation
TStart = tic;
% Before the Disturbances
[time_vector1, x1] = DSA(StoichiometryMatrix, PropensityFunction, NetworkParameters, InitialCondition, t1, N0, Solver);
% First Disturbance at time t1
[time_vector2, x2] = DSA(StoichiometryMatrix, PropensityFunction, NetworkParameters1, x1(:,end), t2-t1, N1, Solver);
% Second Disturbance at time t2
[time_vector3, x3] = DSA(StoichiometryMatrix, PropensityFunction, NetworkParameters2, x2(:,end), tf - t2, N2, Solver);
SimTime_APIF_Deg = toc(TStart);
% Construct Full Solution
time_vector_APIF_Deg = [time_vector1, t1 + time_vector2(2:end), t2 + time_vector3(2:end)];
x_APIF_Deg = [x1, x2(:,2:end), x3(:,2:end)];
clear x1 x2 x3 time_vector1 time_vector2 time_vector3 ...
      NetworkParameters1 NetworkParameters2 NetworkParameters ...
      StoichiometryMatrix PropensityFunction;
  
%% Simulation (APIDF, Degradation)
% Controller Parameters
ControllerParameters = AIFParameters;
ControllerParameters.beta = 50;
ControllerParameters.theta = ControllerParameters.theta + ControllerParameters.beta;
ControllerParameters.delta = 0.1;
ControllerParameters.n = 1;
% Network Parameters
NetworkParameters = CombineStructures(PlantParameters, ControllerParameters);
NetworkParameters1 = CombineStructures(PlantParameters1, ControllerParameters);
NetworkParameters2 = CombineStructures(PlantParameters2, ControllerParameters);
% Initial Condition
InitialCondition = [IC; 0; 0];
% Stoichiometry Matrix and Propensity Function
StoichiometryMatrix = StoichiometryMatrix_APIDF_PositiveZ1_Deg();
PropensityFunction = @PropensityFunction_APIDF_PositiveZ1_Deg;
% Deterministic Simulation
TStart = tic;
% Before the Disturbances
[time_vector1, x1] = DSA(StoichiometryMatrix, PropensityFunction, NetworkParameters, InitialCondition, t1, N0, Solver);
% First Disturbance at time t1
[time_vector2, x2] = DSA(StoichiometryMatrix, PropensityFunction, NetworkParameters1, x1(:,end), t2-t1, N1, Solver);
% Second Disturbance at time t2
[time_vector3, x3] = DSA(StoichiometryMatrix, PropensityFunction, NetworkParameters2, x2(:,end), tf - t2, N2, Solver);
SimTime_APIF_Deg = toc(TStart);
% Construct Full Solution
time_vector_APIDF_Deg = [time_vector1, t1 + time_vector2(2:end), t2 + time_vector3(2:end)];
x_APIDF_Deg = [x1, x2(:,2:end), x3(:,2:end)];
clear x1 x2 x3 time_vector1 time_vector2 time_vector3 ...
      NetworkParameters1 NetworkParameters2 NetworkParameters ...
      StoichiometryMatrix PropensityFunction;

%% Figure Settings
ScalingFactor = 1.3;
Figure_Width = 4 * SS;
Figure_Height = 3.5 * SS;
FontSize = ScalingFactor*5 * SS;
FontSize_Small = ScalingFactor*4 * SS;
LineWidth = ScalingFactor*0.65 * SS;
LineWidth_Thick = ScalingFactor*1 * SS;
MarkerSize = ScalingFactor*3 * SS;
ControllerColor = [152,78,163]/255;
PlantColor = [255,127,0]/255;
DisturbanceColor = [255,215,0]/ 255;
DisturbanceColor2 = [247,129,191] / 255;

%% Set Figure 1
Figure1_Name = 'LDLC';
Handle_Figure1 = figure();
    Handle_Figure1.Color = [1 1 1];
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [-15, 100, 1.5*Figure_Width, Figure_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];
    set(Handle_Figure1, 'InvertHardCopy', 'off'); 
Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Position = [0.15, 0.19, 0.8, 0.75];
    Handle_Axis1.Box = 'on';
    Handle_Axis1.FontSize = FontSize;
    hold(Handle_Axis1, 'on');
    grid(Handle_Axis1, 'on');
    Handle_Axis1.XMinorGrid = 'off';
    Handle_Axis1.YMinorGrid = 'off';
    Handle_Axis1.XLabel.String = 'Time (Days)';
    Handle_Axis1.YLabel.String = 'LDL-C (mg/dL)';
    Handle_Axis1.XScale = 'linear';
    Handle_Axis1.YScale = 'linear'; 
    Handle_Axis1.XLim = [t0, tf];
    Handle_Axis1.XTick = (t1 : 100 : tf);
    Handle_Axis1.XTickLabel = 0 : 100 : tf-t1;
    
%% Set Figure 2
Figure2_Name = 'IC';
Handle_Figure2 = figure();
    Handle_Figure2.Color = [1 1 1];
    Handle_Figure2.PaperUnits = 'centimeters';
    Handle_Figure2.Units = 'centimeters';
    Handle_Figure2.Position = [-15 + Figure_Width, 100, 1.5*Figure_Width, Figure_Height];
    Handle_Figure2.PaperPositionMode = 'auto';
    Handle_Figure2.PaperSize = [Handle_Figure2.PaperPosition(3), Handle_Figure2.PaperPosition(4)];
    set(Handle_Figure2, 'InvertHardCopy', 'off');
Handle_Axis2 = axes(Handle_Figure2);
    Handle_Axis2.Position = [0.15, 0.19, 0.8, 0.75];
    Handle_Axis2.Box = 'on';
    Handle_Axis2.FontSize = FontSize;
    hold(Handle_Axis2, 'on');
    grid(Handle_Axis2, 'on');
    Handle_Axis2.XMinorGrid = 'off';
    Handle_Axis2.YMinorGrid = 'off';
    Handle_Axis2.XLabel.String = 'Time (Days)';
    Handle_Axis2.YLabel.String = 'IC (mg/dL)';
    Handle_Axis2.XScale = 'linear';
    Handle_Axis2.YScale = 'linear'; 
    Handle_Axis2.XLim = [t0, tf];
    Handle_Axis2.XTick = (t1 : 100 : tf);
    Handle_Axis2.XTickLabel = 0 : 100 : tf-t1;
    Handle_Axis2.YLim = [1250, 3400];
    
%% Plotting the Results
% LDL-C
plot(Handle_Axis1, time_vector_OL, x_OL(Output_Index,:), 'LineWidth', LineWidth, 'Color', 'k');
plot(Handle_Axis1, time_vector_AIF, x_AIF(Output_Index,:), 'LineWidth', LineWidth, 'Color', MyBlue);
plot(Handle_Axis1, time_vector_APIF_Deg, x_APIF_Deg(Output_Index,:), 'LineWidth', LineWidth, 'Color', MyRed);
plot(Handle_Axis1, time_vector_APIDF_Deg, x_APIDF_Deg(Output_Index,:), 'LineWidth', LineWidth, 'Color', MyGreen);
plot(Handle_Axis1, [t1 t1], Handle_Axis1.YLim, 'Color', DisturbanceColor, 'LineWidth', LineWidth, 'LineStyle', '--');
plot(Handle_Axis1, [t2 t2], Handle_Axis1.YLim, 'Color', DisturbanceColor2, 'LineWidth', LineWidth, 'LineStyle', '--');
Handle_Legend1 = legend(Handle_Axis1, 'Open Loop', 'I', 'PI', 'PID');
Handle_Legend1.Location = 'North West';
Handle_Legend1.FontSize = FontSize;
% IC
plot(Handle_Axis2, time_vector_OL, x_OL(Input_Index,:), 'LineWidth', LineWidth, 'Color', 'k');
plot(Handle_Axis2, time_vector_AIF, x_AIF(Input_Index,:), 'LineWidth', LineWidth, 'Color', MyBlue);
plot(Handle_Axis2, time_vector_APIF_Deg, x_APIF_Deg(Input_Index,:), 'LineWidth', LineWidth, 'Color', MyRed);
plot(Handle_Axis2, time_vector_APIDF_Deg, x_APIDF_Deg(Input_Index,:), 'LineWidth', LineWidth, 'Color', MyGreen);
plot(Handle_Axis2, [t1 t1], Handle_Axis2.YLim, 'Color', DisturbanceColor, 'LineWidth', LineWidth, 'LineStyle', '--');
plot(Handle_Axis2, [t2 t2], Handle_Axis2.YLim, 'Color', DisturbanceColor2, 'LineWidth', LineWidth, 'LineStyle', '--');
Handle_Legend2 = legend(Handle_Axis2, 'Open Loop', 'I', 'PI', 'PID');
Handle_Legend2.Location = 'East';
Handle_Legend2.FontSize = FontSize;


%% Saving Figures
if Save_Flag == 1  
    Handle_Figure1.Color = 'none';
    Handle_Figure2.Color = 'none';
    Handle_Figure3.Color = 'none';
    Handle_Figure4.Color = 'none';
    print(Handle_Figure1, Figure1_Name,'-dpdf');
    print(Handle_Figure2, Figure2_Name,'-dpdf');
    print(Handle_Figure3, Figure3_Name,'-dpdf');
    print(Handle_Figure4, Figure4_Name,'-dpdf');
end