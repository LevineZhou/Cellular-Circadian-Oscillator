close all; clear all; clc;

% set parameters
% alpha0 = 12;
V_n = 1e-5;
V_c = 2e-4;
r = 5;
beta = 10; % /h
v = (2 * pi) / 22; % /h
gamma_m = v;
gamma_p = v;
delta_m = v;
delta_p = v;
K = 200;

% in the deterministic case
% sin light level
f = @(t, X)[(alpha(t) / V_n) * ((K / (K + X(4))) ^ r) - gamma_m * X(1);
    gamma_m * (V_n / V_c) * X(1) - delta_m * X(2);
    beta * X(2) - gamma_p * X(3);
    gamma_p * (V_c / V_n) * X(3) - delta_p * X(4)];

[t, X] = ode45(f, [0, 720], [10.9 1.88 129.21 3531.75]);

% constant light level
f = @(t1, X1)[(alpha1(t1) / V_n) * ((K / (K + X1(4))) ^ r) - gamma_m * X1(1);
    gamma_m * (V_n / V_c) * X1(1) - delta_m * X1(2);
    beta * X1(2) - gamma_p * X1(3);
    gamma_p * (V_c / V_n) * X1(3) - delta_p * X1(4)];

[t1, X1] = ode45(f, [0, 720], [10.9 1.88 129.21 3531.75]);

% constant light level
f = @(t4, X4)[(alpha4(t4) / V_n) * ((K / (K + X4(4))) ^ r) - gamma_m * X4(1);
    gamma_m * (V_n / V_c) * X4(1) - delta_m * X4(2);
    beta * X4(2) - gamma_p * X4(3);
    gamma_p * (V_c / V_n) * X4(3) - delta_p * X4(4)];

[t4, X4] = ode45(f, [0, 720], [10.9 1.88 129.21 3531.75]);

% all day light
f = @(t2, X2)[(alpha2(t2) / V_n) * ((K / (K + X2(4))) ^ r) - gamma_m * X2(1);
    gamma_m * (V_n / V_c) * X2(1) - delta_m * X2(2);
    beta * X2(2) - gamma_p * X2(3);
    gamma_p * (V_c / V_n) * X2(3) - delta_p * X2(4)];

[t2, X2] = ode45(f, [0, 720], [10.9 1.88 129.21 3531.75]);

% change time zone
f = @(t3, X3)[(alpha3(t3) / V_n) * ((K / (K + X3(4))) ^ r) - gamma_m * X3(1);
    gamma_m * (V_n / V_c) * X3(1) - delta_m * X3(2);
    beta * X3(2) - gamma_p * X3(3);
    gamma_p * (V_c / V_n) * X3(3) - delta_p * X3(4)];

[t3, X3] = ode45(f, [0, 1440], [10.9 1.88 129.21 3531.75]);

f = @(t5, X5)[(alpha5(t5) / V_n) * ((K / (K + X5(4))) ^ r) - gamma_m * X5(1);
    gamma_m * (V_n / V_c) * X5(1) - delta_m * X5(2);
    beta * X5(2) - gamma_p * X5(3);
    gamma_p * (V_c / V_n) * X5(3) - delta_p * X5(4)];

[t5, X5] = ode45(f, [0, 1440], [10.9 1.88 129.21 3531.75]);

% plot the quantities to see the oscillation
% figure(1);
% plot(t1, X1(:, 1));
% hold on
% plot(t1, X1(:, 2));
% hold on
% plot(t1, X1(:, 3));
% hold on
% plot(t1, X1(:, 4));
% grid on
% legend('M_n','M_c', 'P_c', 'P_n', 'Location', 'best')
% xlabel('time(h)'); ylabel('/pL')

test_t = [0, 1440];

light_level0 = zeros(1440, 1);
light_time0 = zeros(1440, 1);
for i = 1:1440
    light_level0(i) = alpha5(i);
    light_time0(i) = i;
end


light_level = zeros(1440, 1);
light_time = zeros(1440, 1);
for i = 1:1440
    light_level(i) = alpha3(i);
    light_time(i) = i;
end

subplot(2, 1, 1)
plot(t1, X1(:, 1));
hold on
plot(t1, X1(:, 2));
hold on
plot(t1, X1(:, 3));
hold on
plot(t1, X1(:, 4));
grid on
legend('M_n','M_c', 'Location', 'best')
xlabel('time(h)'); ylabel('/pL')
xlim([600, 720])
title('mRNA')

% subplot(2, 1, 2)
% plot(t1, X1(:, 1));
% hold on
% plot(t1, X1(:, 2));
% hold on
% plot(t1, X1(:, 3));
% hold on
% plot(t1, X1(:, 4));
% grid on
% legend('P_c', 'P_n', 'Location', 'best')
% xlabel('time(h)'); ylabel('/pL')
% xlim([600, 720])
% title('protein')

% subplot(2, 1, 2)
% plot(light_time, light_level);
% legend('light', 'Location', 'best')
% xlabel('time(h)'); ylabel('light')
% ylim([0, 30])
% title('light level')

% figure(2);
% plot(t, X(:, 1));
% hold on
% plot(t, X(:, 2));
% grid on
% legend('M_n','M_c', 'Location', 'best')
% xlabel('time(h)'); ylabel('/pL')
% 
% figure(3);
% plot(t, X(:, 3));
% hold on
% plot(t, X(:, 4));
% grid on
% legend('P_c', 'P_n', 'Location', 'best')
% xlabel('time(h)'); ylabel('/pL')

% figure(4);
% plot(t, X(:, 3));
% hold on
% plot(t1, X1(:, 3));
% hold on
% plot(t2, X2(:, 3));
% grid on
% legend('continuous dark', '12:12', 'constant light', 'Location', 'best')
% xlim([600, 720])
% xlabel('time(h)'); ylabel('Pc/pL')

% figure(5);
% plot(t1, X1(:, 3));
% hold on
% plot(t3, X3(:, 3));
% grid on
% legend('12:12', 'jet lag 12', 'Location', 'best')
% xlabel('time(h)'); ylabel('Pc/pL')
% xlim(test_t)

% subplot(2, 2, 1)
% plot(t1, X1(:, 3));
% legend('12:12', 'Location', 'best')
% xlabel('time(h)'); ylabel('Pc/pL')
% xlim(test_t)
% title('Origin')
% 
% subplot(2, 2, 2)
% plot(t4, X4(:, 3));
% legend('12:12', 'Location', 'best')
% xlabel('time(h)'); ylabel('Pc/pL')
% xlim(test_t)
% title('Destination')

% subplot(2, 2, 1)
% plot(t5, X5(:, 3));
% legend('jet lag', 'Location', 'best')
% xlabel('time(h)'); ylabel('Pc/pL')
% xlim([600, 1080]);
% title('the cell LN')
% 
% subplot(2, 2, 2)
% plot(t3, X3(:, 3));
% legend('jet lag', 'Location', 'best')
% xlabel('time(h)'); ylabel('Pc/pL')
% xlim([600, 1080]);
% title('the cell NL')
% 
% subplot(2, 2, 3)
% plot(light_time0, light_level0);
% legend('light', 'Location', 'best')
% xlabel('time(h)'); ylabel('light')
% xlim([600, 1080]); 
% ylim([0, 30])
% title('light level')
% 
% subplot(2, 2, 4)
% plot(light_time, light_level);
% legend('light', 'Location', 'best')
% xlabel('time(h)'); ylabel('light')
% xlim([600, 1080]); 
% ylim([0, 30])
% title('light level')

function f = alpha(t)
    alpha_0 = 12;
    alpha_bar = 6;
	f = alpha_0 + alpha_bar * (1 + epsilon(t));
%     f = alpha_0;
end

function e = epsilon(t)
	e = sin((2*pi/24)*t);
end

function f = alpha1(t)
    alpha_0 = 12;
    alpha_bar = 6;
    if (mod(t, 24) <= 12)
        f = alpha_0 + alpha_bar * 2;
    else
        f = alpha_0;
    end
end

function f = alpha4(t)
    alpha_0 = 12;
    alpha_bar = 6;
    if (mod(t, 24) <= 12)
        f = alpha_0;
    else
        f = alpha_0 + alpha_bar * 2;
    end
end

function f = alpha2(t)
    alpha_0 = 12;
    alpha_bar = 6;
    if (t >= 0)
        f = alpha_0 + alpha_bar * 2;
    end
end

% jet lag simualtion
function f = alpha3(t)
    alpha_0 = 12;
    alpha_bar = 6;
    T = 716;
    % stay in the current time zone
    if (t < T)
        if (mod(t, 24) <= 12)
            f = alpha_0;
        else
            f = alpha_0 + alpha_bar * 2;
        end
%         f = alpha(t);
    end
    % start travelling
    if (t <= T + 7) && (t >= T)
        f = alpha_0;
    end
    % end traveling and in a new time zone
    if (t > T + 7)
        if (mod(t + 5, 24) <= 12)
            f = alpha_0 + alpha_bar * 2;
        else
            f = alpha_0;
        end
%         f = alpha(t + 5);
    end
end

function f = alpha5(t)
    alpha_0 = 12;
    alpha_bar = 6;
    T = 711;
    % stay in the current time zone
    if (t < T)
        if (mod(t + 5, 24) <= 12)
            f = alpha_0;
        else
            f = alpha_0 + alpha_bar * 2;
        end
%         f = alpha(t + 5);
    end
    % start travelling
    if (t <= T + 9) && (t >= T)
        f = alpha_0;
    end
    % end traveling and in a new time zone
    if (t > T + 9)
        if (mod(t, 24) <= 12)
            f = alpha_0 + alpha_bar * 2;
        else
            f = alpha_0;
        end
%         f = alpha(t);
    end
end
