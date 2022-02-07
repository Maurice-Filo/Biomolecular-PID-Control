function Cost = Compute_PerformanceIndex(t, x, r, Threshold, Weights)
% Weights = [Settling Time, Rise Time, Overshoot]
Cost = Weights(1) * Compute_SettlingTime(t, x, r, Threshold) + Weights(2) * Compute_RiseTime(t, x, r, Threshold) + Weights(3) * Compute_Overshoot(x, r);
end

