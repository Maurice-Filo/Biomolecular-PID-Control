function T = Compute_RiseTime(t, x, r, Threshold)
N = length(x);
T = -1;
for i = 2 : N
    if (x(i) > r * (1 - Threshold)) && (x(i-1) <= r * (1 - Threshold))
        T = t(i);
        break;
    end
end
end

