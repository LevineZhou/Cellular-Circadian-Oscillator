close all; clear all; clc;

f = @(t, X)[(alpha(t) / V_n) * ((K / (K + X(4))) ^ r) - gamma_m * X(1);
    gamma_m * (V_n / V_c) * X(1) - delta_m * X(2);
    beta * X(2) - gamma_p * X(3);
    gamma_p * (V_c / V_n) * X(3) - delta_p * X(4)];