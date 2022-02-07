function T = Compute_SettlingTime(t, x, r, Threshold)
N = length(x);
T = -1;
for i = N : -1 : 1
    if (x(i) <= r - Threshold*r) || (x(i) >= r + Threshold*r)
        T = t(i);
        break;
    end
end
end

