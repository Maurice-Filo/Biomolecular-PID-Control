function [C, t_vector] = ComputeCovar(X1, X2, T, N)
M = size(T);
for i = 1 : M
    X{i} = X1{i} .* X2{i};
end
[t_vector, MeanX, ~] = ComputeStatistics(T, X, N);
[~, MeanX1, ~] = ComputeStatistics(T, X1, N);
[~, MeanX2, ~] = ComputeStatistics(T, X2, N);
C = MeanX - MeanX1 .* MeanX2;
end

