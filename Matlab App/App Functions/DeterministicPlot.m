function DeterministicPlot(time_vector, x, Parameters, Title, Names)
    %% Create Colors for Plant and Controller Species
    Colors = lines(2);
    %% Total Number of Species
    N = size(x,1);
    %% Create Figure and Subplots 
    figure();
    Handle = gobjects(N,1);
    Nr = round(sqrt(N));
    if round(sqrt(N))^2 >= N
        Nc = Nr;
    else 
        Nc = Nr+1;
    end
    for k = 1 : N
        Handle(k) = subplot(Nr, Nc, k);
    end
    %% Main Title
    sgtitle([Title, ', Deterministic Setting']);
    %% Plotting...
    for k = 1 : N
        if N-k < 2
            c = Colors(1,:);
        else 
            c = Colors(2,:);
        end
        plot(Handle(k), time_vector, x(k,:), 'Color', c, 'LineWidth', Parameters.LineWidth);
        Handle(k).XLabel.String = '$t$';
        Handle(k).XLabel.Interpreter = 'latex';
        Handle(k).YLabel.String = ['$', Names{k}, '$(t)'];
        Handle(k).YLabel.Interpreter = 'latex';
        Handle(k).XLim = [0 time_vector(end)];
        Handle(k).XGrid = 'on';
        Handle(k).YGrid = 'on';
        Handle(k).XMinorGrid = 'on';
        Handle(k).YMinorGrid = 'on';
        Handle(k).FontSize = Parameters.FontSize;
        Handle(k).FontName = 'Times New Roman';
    end
end

