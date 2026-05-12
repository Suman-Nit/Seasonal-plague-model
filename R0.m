
close all; clc;

tspan = linspace(0, 365, 1000);
r_R = 5/365; p = 0.65; K_R = 250000; d_R = 0.00055;  m_R = 0.0547; g_R = 0.06; mu_R = 0.03; m_H = 0.1; g_H = 0.34;
a = 3/250000; r_F = 0.0084; mu_F = 0.008; r_H = 0.00011; d_H = 0.00011; N_H=10000;
sigma1 = 0.150994; sigma2 = 0.033235; k_f=0.019684; alpha =  0.64195; omega_1 = 160;
beta_R_bar = 0.53057; beta_H_bar = 0.0097422;  K_F_bar = 3; d_F_bar = 0.20569;


target_eigenvalue = 1;

lambda_initial_guess = 2;


options = optimset('TolFun', 1e-10, 'TolX', 1e-10);

lambda_solution = fzero(@(lambda_guess) model_lambda_time_avg(lambda_guess, tspan, beta_R_bar, sigma1, beta_H_bar, sigma2,...
    d_F_bar, alpha, K_F_bar, k_f, d_R, m_R, a, K_R, g_R, N_H, d_H, m_H, omega_1) - target_eigenvalue, ...
    lambda_initial_guess, options);

lambda_values = lambda_solution
