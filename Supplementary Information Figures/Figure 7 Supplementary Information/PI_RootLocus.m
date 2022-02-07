%% Clear Workspace
close all
clc
clear;
Save_Flag = 0;

%% Parameters
gamma_1 = 0.1;
gamma_2 = 0.1;
k_1 = 0.1;
K_S = 2;
N = 5;
K_P_vector = linspace(0, 50, 1000);
K_I_Threshold = gamma_1 * gamma_2 * (gamma_1 + gamma_2) / k_1 / K_S;

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
PIColor1 = [0, 128, 128]/255; % Green
PIColor2 = [230, 25, 75]/255; % Red
PIColor3 = [0, 0, 128]/255; % Blue

%% Set Figure 1
Figure1_Name = 'PI_RootLocus';
Handle_Figure1 = figure();
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [-30, 10 + 3.5*Figure_Height, Figure_Width, Figure_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];

%% Axis for Root Locus Starting from Stable Poles
Handle_Axis1 = subplot(1,2,1);
    Handle_Axis1.Position = [0.07, 0.2, 0.4, 0.75];
    hold(Handle_Axis1, 'on');
K_I_vector1 = linspace(0, K_I_Threshold, N);
for i = 1 : length(K_I_vector1)
    K_I = K_I_vector1(i);
    System = tf([k_1 0], [1, gamma_1+gamma_2, gamma_1*gamma_2, k_1*K_S*K_I]);
    Handle = rlocusplot(Handle_Axis1, System, K_P_vector);
    Options = getoptions(Handle);
        Options.Title.String = '';
        Options.Title.Interpreter = 'latex';
        Options.Title.FontSize = FontSize;
        Options.XLabel.String = '';
        Options.XLabel.FontSize = FontSize/10000;
        Options.YLabel.String = '';
        Options.YLabel.FontSize = FontSize/10000;
        setoptions(Handle, Options);
end
Handle_Axis1.XColor = 'k';
Handle_Axis1.YColor = 'k';
Handle_Axis1.Box = 'on';
Handle_Axis1.BoxStyle = 'full';
Handle_Axis1.LineWidth = LineWidth_Thin;
Handle_Axis1.FontSize = FontSize;
Handle_Axis1.XLim = [-0.7, 0.3];
Handle_Axis1.YLim = [-2.2, 2.2];
Handle_Axis1.XLabel.String = 'Real Axis';
Handle_Axis1.YLabel.String = 'Imaginary Axis';
Handle_Axis1.YLabel.Color = 'k';
Handle_Axis1.XLabel.Color = 'k';
Handle_Axis1.Title.String = '$(\gamma_1 + \gamma_2) \gamma_1 \gamma_2 > k_1 K_S K_I$';
Handle_Axis1.Title.Interpreter = 'latex';
Handle_Axis1.Title.Color = 'k';
Handle_Axis1.XTick = -0.6 : 0.3 : 0.3;
Handle_Annotation1 = annotation(Handle_Figure1, 'textarrow', [0.214331597222222,0.314331597222222], [0.352916666666667,0.202916666666667]);
    Handle_Annotation1.String = '$-\frac{\gamma_1 + \gamma_2}{2}$';
    Handle_Annotation1.Interpreter = 'latex';
    Handle_Annotation1.FontSize = FontSize*1.25;
   

%% Axis for Root Locus Starting from Unstable Poles
Handle_Axis2 = subplot(1,2,2);
    Handle_Axis2.Position = [0.57, 0.2, 0.4, 0.75];
    hold(Handle_Axis2, 'on');
K_I_vector2 = linspace(K_I_Threshold, 100 * K_I_Threshold, N);
for i = 1 : length(K_I_vector2)
    K_I = K_I_vector2(i);
    System = tf([k_1 0], [1, gamma_1+gamma_2, gamma_1*gamma_2, k_1*K_S*K_I]);
    Handle = rlocusplot(Handle_Axis2, System, K_P_vector);
    Options = getoptions(Handle);
    if i == length(K_I_vector1)
        Options.Title.String = '';
        Options.Title.Interpreter = 'latex';
        Options.Title.FontSize = FontSize;
        Options.XLabel.String = '';
        Options.XLabel.FontSize = FontSize/10000;
        Options.YLabel.String = '';
        Options.YLabel.FontSize = FontSize/10000;
        setoptions(Handle, Options);
    end
end
Handle_Axis2.XColor = 'k';
Handle_Axis2.YColor = 'k';
Handle_Axis2.Box = 'on';
Handle_Axis2.BoxStyle = 'full';
Handle_Axis2.LineWidth = LineWidth_Thin;
Handle_Axis2.FontSize = FontSize;
Handle_Axis2.XLim = [-0.7, 0.3];
Handle_Axis2.YLim = [-2.2, 2.2];
Handle_Axis2.XLabel.String = 'Real Axis';
Handle_Axis2.YLabel.String = 'Imaginary Axis';
Handle_Axis2.YLabel.Color = 'k';
Handle_Axis2.XLabel.Color = 'k';
Handle_Axis2.Title.String = '$(\gamma_1 + \gamma_2) \gamma_1 \gamma_2 < k_1 K_S K_I$';
Handle_Axis2.Title.Interpreter = 'latex';
Handle_Axis2.Title.Color = 'k';
Handle_Axis2.XTick = -0.6 : 0.3 : 0.3;
Handle_Annotation2 = annotation(Handle_Figure1, 'textarrow', [0.716076388888886,0.816076388888887], [0.349976851851852,0.199976851851853]);
    Handle_Annotation2.String = '$-\frac{\gamma_1 + \gamma_2}{2}$';
    Handle_Annotation2.Interpreter = 'latex';
    Handle_Annotation2.FontSize = FontSize*1.25;
    
%% Save Figures
if Save_Flag == 1
    Handle_Figure1.Color = 'none';
    set(Handle_Figure1, 'InvertHardCopy',  'off');
    print(Handle_Figure1, Figure1_Name, '-dpdf', '-painters');
end