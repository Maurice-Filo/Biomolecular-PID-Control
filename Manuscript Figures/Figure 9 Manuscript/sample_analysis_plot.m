clear all;
close all;
clc;

%% Parameters
exptIndx{1} = 'I_168cells'; 
exptIndx{2} = 'PI_128cells'; 
exptIndx{3} = 'PID_178cells'; 
setPoint = 14;
numFrames = 121;
startRelevantFrameFFT = 30;
endRelevantFrameFFT = 121;
samplingTime = 2*60; % seconds
fullTimeAxis = [0:(numFrames-1)]*samplingTime;
colorsExp = [
     55,126,184
     228,56,58
     77,175,74
    ]/255;
plotsVertical = 2;
plotsHorizontal = 1;

numExpt = numel(exptIndx);
relevantIndxCount = zeros(numExpt,1);

% Loop through experimental data
for loopNum=1:numExpt
    
    load(['source_data' filesep 'data_' exptIndx{loopNum} '.mat']);
    [relevantIndxCount(loopNum),~] = size(nascentRNAcount);        
    
    %% Compute FFT (look help fft!)
    zeroMeanOutput = nascentRNAcount - setPoint;    
    lengthSignal = endRelevantFrameFFT - startRelevantFrameFFT + 1;
    samplingFrequency = 1/(samplingTime);
    P_mat = [];
    for cellItr = 1:relevantIndxCount(loopNum)
        [freqAxis,P1] = fft_cells(zeroMeanOutput(cellItr,startRelevantFrameFFT:endRelevantFrameFFT),lengthSignal,samplingFrequency);
        P_mat = [P_mat;P1];
    end
    meanP = mean(P_mat);
    
    %% Plot FFT
    figure(1);

    psd_stdError = [];
    for i=1:size(P_mat,2)
        psd_stdError = [psd_stdError, std(P_mat(:,i))/sqrt(relevantIndxCount(loopNum))];
    end
    psd_errbar = [psd_stdError ; psd_stdError];

    subplot(plotsVertical,plotsHorizontal,2);
    hold on;
    psd_s1(loopNum) = shadedErrorBar(freqAxis(2:end).*(2*pi*samplingTime), meanP(2:end), psd_errbar(:,2:end));
    psd_s1(loopNum).mainLine.LineWidth = 2;
    psd_s1(loopNum).mainLine.Color = colorsExp(loopNum,:);
    psd_s1(loopNum).patch.FaceColor = colorsExp(loopNum,:);    
    hold off; 
        
    %% Plot Mean
    meanNascentRNAcount = mean(nascentRNAcount);
    stdError = [];
    for i=1:numFrames
        stdError = [stdError, std(nascentRNAcount(:,i))/sqrt(relevantIndxCount(loopNum))];
    end
    errbar = [stdError ; stdError];
    
    subplot(plotsVertical,plotsHorizontal,1);    
    hold on;
    s1(loopNum) = shadedErrorBar(fullTimeAxis./60, meanNascentRNAcount, errbar); % converting time axis from seconds to minutes
    s1(loopNum).mainLine.LineWidth = 2;
    s1(loopNum).mainLine.Color = colorsExp(loopNum,:);
    s1(loopNum).patch.FaceColor = colorsExp(loopNum,:);    
    s2(loopNum) = scatter(fullTimeAxis./60, meanNascentRNAcount, 'MarkerEdgeColor',colorsExp(loopNum,:),...
              'MarkerFaceColor',colorsExp(loopNum,:),...
              'LineWidth',1);
    s2(loopNum).SizeData = 10;
    hold off; 
    
end

figure(1);
subplot(plotsVertical,plotsHorizontal,2);
title('Mean of the PSDs of single cell output trajectories')
xlabel('\omega (rad./samples)')
ylabel('S(\omega)')
legend([psd_s1(1).mainLine psd_s1(2).mainLine psd_s1(3).mainLine],{'I','PI','PID'})
set(gca, 'FontSize', 13)
xlim([0 pi/2])
grid on;

figure(1);
subplot(plotsVertical,plotsHorizontal,1); hold on;
title('Mean of single cell output trajectories')
xlabel('Time (minutes)')
ylabel('Nascent RNA count (in vivo)')

ref = refline(0,setPoint);
ref.LineStyle = '--';
ref.LineWidth = 2;
ref.Color = 'k';
legend([s1(1).mainLine s1(2).mainLine s1(3).mainLine ref],{'I','PI','PID','Set point'})
set(gca, 'FontSize', 13)
xlim([0 240])
hold off;
grid on;


function [f,P1] = fft_cells(x,L,Fs)
n = 2^(nextpow2(L)+3);
Y = fft(x,n);
P2 = abs(Y/n).^2;
P1 = P2(1:n/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(n/2))/n;
end
