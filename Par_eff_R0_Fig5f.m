
clear; close all; clc;


tspan = linspace(0,365,1000);
K_R = 250000; a = 3/K_R; d_R = 0.00055; m_R = 0.0547; g_R = 0.06;
N_H = 10000; d_H = 0.00011; m_H = 0.1; omega_1 = 160;
sigma1 = 0.150994; beta_R_bar = 0.53057;
 sigma2 = 0.033235;  
  beta_H_bar = 0.0097422;  
% d_F_bar = 0.20569;
 alpha =  0.64195; 
  k_f=0.019684; % K_F_bar = 3;


% k_f_vals = linspace(0,1,20);
% K_F_bar_vals = linspace(1,15,20);
K_F_bar_vals = linspace(2,11,20);
 d_F_bar_vals = linspace(0,1,20);
[BR, S1] = meshgrid(d_F_bar_vals,K_F_bar_vals);
R0 = nan(size(BR));

opts = optimset('TolFun',1e-8,'TolX',1e-8,'Display','off');
target = 1; init_guess = 2 ;

for i = 1:size(BR,1)
    last_guess = init_guess;
    for j = 1:size(BR,2)
        br = BR(i,j);
        s1 = S1(i,j);
        f = @(lam) model_lambda(lam, tspan, beta_R_bar, sigma1, beta_H_bar, sigma2, ...
            br, alpha, s1, k_f, d_R, m_R, a, K_R, g_R, N_H, d_H, m_H, omega_1) - target;
            sol = fzero(f, last_guess, opts);
        R0(i,j) = sol;
    end
end

figure;
contourf(BR, S1, R0, 20, 'LineColor','none'); colorbar; hold on;
plot(0.20569, 3, 'v', 'MarkerSize',10, 'MarkerFaceColor','k', 'MarkerEdgeColor','k');
xlabel('$\overline{d_F}$','Interpreter','latex'); ylabel('$\overline{K_F}$','Interpreter','latex');
 title('R_0 contour');
 colormap(jet);
