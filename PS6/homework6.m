clear;clc;

global rho sigma beta r

rho = 0.92;
sigma_e = 0.2;
beta = 0.95;
r = 0.025;

num_st = 7     % number of state
sigma_y = sigma_e / sqrt(1 - rho^2);

% construct the transition matrix and state vector
[P_matrix, y] = rouwen(rho, 0, sigma_y, num_st);

num_sm = 1000; % number of simulated periods

ydiscR = hitm_z(y, P_matrix,rand(num_sm,1));
plot(ydiscR)
xlabel('Period, t')
ylabel('ln(y_t)')
title('Rowenhorst, \rho=0.92 and \sigma \epsilon =0.2')
T=autocorr(ydiscR,1);
V=std(ydiscR);
