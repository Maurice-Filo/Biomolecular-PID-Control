function StochasticPlot(T, X, Parameters, Title, Names)
    %% Create Colors for Plant and Controller Species
    Colors = lines(2);
    %% Total Number of Species
    N = size(X{1},1);
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
    sgtitle([Title, ', Stochastic Setting']);
    %% Plotting...
    for k = 1 : N
        if N-k < 2
            c = Colors(1,:);
        else 
            c = Colors(2,:);
        end
        hold(Handle(k), 'on');
        for i = 1 : length(T)
            plot(Handle(k), T{i}, X{i}(k,:), 'Color', [c, 0.2], 'LineWidth', Parameters.LineWidth);
        end
        Handle(k).XLabel.String = '$t$';
        Handle(k).XLabel.Interpreter = 'latex';
        Handle(k).YLabel.String = ['$' Names{k}, '$(t)'];
        Handle(k).YLabel.Interpreter = 'latex';
        Handle(k).XGrid = 'on';
        Handle(k).YGrid = 'on';
        Handle(k).XMinorGrid = 'on';
        Handle(k).YMinorGrid = 'on';
        Handle(k).FontSize = Parameters.FontSize;
        Handle(k).FontName = 'Times New Roman';
    end
end

