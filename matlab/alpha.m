close all; clear all; clc;

% sin function for light level
function f = alpha(t)
	alpha_0 = 0.16;
    alpha_bar = 0.02;
	f = alpha_0 + alpha_bar * (1 + epsilon(t));
end

function e = epsilon(t)
	e = sin((2*pi/24)*t);
end

% 12:12 light-dark cycles
function f = alpha(t)
	alpha_0 = 0.16;
	alpha_bar = 0.02;

	if (t <= 12)
		f = alpha_0 + alpha_bar * 2;
	else
		f = alpha_0 + alpha_bar * 1;
	end
end