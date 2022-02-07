function H = APIDF2_Deg_TF(Gains, P)
s = tf([1, 0], 1);
F = Gains.omega_c / (s + Gains.omega_c);
H = (Gains.K_F + Gains.K_I/s) * P * F / (1 + (Gains.K_P + Gains.K_I*Gains.K_S/s + Gains.K_D*s) * F * P);
H = minreal(H);
end

