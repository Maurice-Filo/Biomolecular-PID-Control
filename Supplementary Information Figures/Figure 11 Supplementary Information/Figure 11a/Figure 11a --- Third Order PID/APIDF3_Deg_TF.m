function H = APIDF3_Deg_TF(Gains, P)
s = tf([1, 0], 1);
F = Gains.omega_0 / (s + Gains.omega_0);
H = (Gains.K_I/s) * P / (1 + ((Gains.K_P + Gains.K_D*s) * F + Gains.K_I*Gains.K_S/s) * P);
H = minreal(H);
end

