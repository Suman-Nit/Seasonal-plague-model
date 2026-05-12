
close all; clc; clear;

r_R = 5/365; p = 0.65; K_R = 250000; d_R = 0.00055; 
m_R = 0.0547; g_R = 0.06; mu_R = 0.03; a = 3/250000; r_F = 0.0084; 
 mu_F = 0.008; r_H = 0.00011; d_H = 0.00011; 
m_H = 0.1; g_H = 0.34;  beta_R_bar = 0.53057;
 N_H=10000; omega_1 = 160;  sigma1 = 0.150994; K_F_bar = 3;
 sigma2 = 0.033235;  alpha =  0.64195;  % k_f=0.019684;
 beta_H_bar = 0.0097422;   d_F_bar = 0.20569; 
 
omega = 365; t_final = 365;

num_s_span = 5000;
s_span_full = linspace(0, 3*omega, num_s_span);
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);


num_beta = 50;                      
k_f_vec = linspace(0, 1, num_beta);  
tau_res = 500;                       
tau = linspace(0, omega, tau_res);

% Storage: rows -> beta, cols -> tau
P_outbreak_mat1 = zeros(num_beta, tau_res);
P_outbreak_mat2 = zeros(num_beta, tau_res);

for j = 1:num_beta
    k_f =k_f_vec(j);

    % time-varying parameter functions 
    beta_R = @(t) beta_R_bar * (1 + sigma1 * cos(2 * pi * t / 365));
    beta_H = @(t) beta_H_bar * (1 + sigma2 * cos(2 * pi * t / 365));
    d_F = @(t) d_F_bar * (1 + alpha * cos(2 * pi * (t  - omega_1)/ 365));
    K_F1 = @(t) K_F_bar * (1 + k_f) .^ (sin(2 * pi * t / 365));

    g1 = @(x1, x2, x3, tau) ((d_R+m_R) + (d_R+m_R*(1-g_R))*K_F1(tau).*x1.*x2) ./ ...
                            ((d_R+m_R) + (d_R+m_R*(1-g_R))*K_F1(tau));
    g2 = @(x1, x2, x3, tau) (beta_R(tau).*(1-exp(-a*K_R)).*x1.*x2 + d_F(tau) + ...
           (1-exp(-a*K_R)) + beta_H(tau).*N_H.*exp(-a*K_R).*x2.*x3 ) ./ ...
           (beta_R(tau).*(1-exp(-a*K_R)) + d_F(tau) + ...
           (1-exp(-a*K_R)) + beta_H(tau).*N_H.*exp(-a*K_R));
    g3 = @(x1, x2, x3, tau) 1;

    
    dHdt = @ (s, H, t) [ - ((d_R+m_R) + (d_R+m_R*(1-g_R))*K_F1(t-s)).*...
               (H(1) - g1(H(1), H(2), H(3), t - s));...
             - (beta_R(t-s).*(1-exp(-a*K_R)) + d_F(t-s) + (1-exp(-a*K_R)) + beta_H(t-s).*N_H.*exp(-a*K_R)).*...
             (H(2) - g2(H(1), H(2), H(3), t - s));...
             -(d_H+m_H).*(H(3) -  g3(H(1), H(2), H(3), t - s))];

    H0 = [0;0;0];
    [~, H] = ode45(@(s,H) dHdt(s, H, t_final), s_span_full, H0, options);

    tol = 1e-6;
    idx_periodic = s_span_full >= 2*omega - tol & s_span_full <= 3*omega + tol;
    s_periodic = s_span_full(idx_periodic) - 2*omega;
    H_periodic = H(idx_periodic, :);

    H1_periodic = @(s) interp1(s_periodic, H_periodic(:,1), mod(s, omega), 'pchip', 'extrap');
    H2_periodic = @(s) interp1(s_periodic, H_periodic(:,2), mod(s, omega), 'pchip', 'extrap');

    for k = 1:tau_res
        ttt = tau(k);
        P_ext_val1 = H1_periodic(omega - ttt);
        P_ext_val2 = H2_periodic(omega - ttt);
        P_outbreak_mat1(j, k) = 1 - P_ext_val1;
        P_outbreak_mat2(j, k) = 1 - P_ext_val2;
    end
end


[T, B] = meshgrid(tau, k_f_vec);   
figure(1);
[c,h]=contourf(T, B, P_outbreak_mat1, 20, 'LineColor','none','ShowText','on');
clabel(c,h, 'FontSize', 10)
pcolor (P_outbreak_mat1); shading interp;
colorbar; colormap jet; drawnow
xlabel('t_0 (day)');
ylabel('$k_f$','Interpreter','latex','FontSize',14);

figure(2);
[c,h]=contourf(T, B, P_outbreak_mat2, 20, 'LineColor','none','ShowText','on');
colorbar;
clabel(c,h, 'FontSize', 10)
pcolor (P_outbreak_mat2); shading interp;
colorbar; colormap jet; drawnow
 xlabel('t_0 (day)'); 
ylabel('$k_f$','Interpreter','latex','FontSize',14);

