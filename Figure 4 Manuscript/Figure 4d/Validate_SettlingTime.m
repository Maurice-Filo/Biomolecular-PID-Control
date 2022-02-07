function Validate_SettlingTime(t, x, r, Threshold)
T = Compute_SettlingTime(t, x, r, Threshold);
figure();
Handle_Axis1 = subplot(1,1,1);
hold on;
plot(Handle_Axis1, t, x, 'Color', 'b', 'LineWidth', 2);
plot(Handle_Axis1, [t(1), t(end)], r*(1 + Threshold)*[1, 1], 'LineWidth', 1, 'Color', 'r');
plot(Handle_Axis1, [t(1), t(end)], r*(1 - Threshold)*[1, 1], 'LineWidth', 1, 'Color', 'r');
plot(Handle_Axis1, T*[1, 1], Handle_Axis1.YLim, 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');
Handle_Axis1.XTick = [0, T, t(end)];
end

