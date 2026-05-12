
clear; close all; clc;

r_R = 5/365; p = 0.65; K_R = 250000; d_R = 0.00055;
m_R = 0.0547; g_R = 0.06; mu_R = 0.03; a = 3/250000; r_F = 0.0084;
mu_F = 0.008; r_H = 0.00011; d_H = 0.00011;
m_H = 0.1; g_H = 0.34; omega_1 = 160; N_H = 10000;

sigma1 = 0.150994; k_f=0.019684;  sigma2 = 0.033235; 
 beta_H_bar = 0.0097422;  beta_R_bar = 0.53057;
  d_F_bar = 0.20569;  K_F_bar = 3; 

omega = 365;

 alpha_values = linspace(0, 1, 50);       
tau_grid  = linspace(0, omega, 500);  

Pmat = zeros(numel(alpha_values), numel(tau_grid));

num_points = 8000;                    
s_span = linspace(0, 3*omega, num_points);
opts = odeset('RelTol',1e-7,'AbsTol',1e-9);

for ii = 1:numel( alpha_values)
     alpha =  alpha_values(ii);
    
    beta_R = @(t) beta_R_bar * (1 + sigma1 * cos(2*pi*t/omega));
    beta_H = @(t) beta_H_bar * (1 + sigma2 * cos(2*pi*t/omega));
    d_F    = @(t) d_F_bar * (1 + alpha * cos(2*pi*(t - omega_1)/omega));
    K_F1   = @(t) K_F_bar * (1 + k_f) .^ (sin(2*pi*t/omega));   
    
    
    g1 = @(x1,x2,x3,tt) ((d_R+m_R) + (d_R+m_R*(1-g_R)).*K_F1(tt).*x1.*x2) ./ ...
                       ((d_R+m_R) + (d_R+m_R*(1-g_R)).*K_F1(tt));
    g2 = @(x1,x2,x3,tt) ( beta_R(tt).*(1-exp(-a*K_R)).*x1.*x2 + d_F(tt) + (1-exp(-a*K_R)) ...
                       ) ./ ...
                       ( beta_R(tt).*(1-exp(-a*K_R)) + d_F(tt) + (1-exp(-a*K_R)) ...
                       + beta_H(tt).*N_H.*exp(-a*K_R) );
    
   
    dHdt = @(s, H, tabs) [ - ((d_R+m_R) + (d_R+m_R*(1-g_R))*K_F1(tabs - s)) .* ...
                              (H(1) - g1(H(1), H(2), 0, tabs - s)); ...
                          - ( beta_R(tabs - s).*(1-exp(-a*K_R)) + d_F(tabs - s) + (1-exp(-a*K_R)) ...
                              + beta_H(tabs - s).*N_H.*exp(-a*K_R) ) .* ...
                              (H(2) - g2(H(1), H(2), 0, tabs - s)) ];
    
    tabs = omega;
    H0 = [0; 0];
    [~, Hsol] = ode45(@(s,H) dHdt(s,H,tabs), s_span, H0, opts);
    
    tol = 1e-9;
    idxP = s_span >= 2*omega - tol & s_span <= 3*omega + tol;
    sP = s_span(idxP) - 2*omega;        
    HP = Hsol(idxP,:);                 
    
    taus_to_eval = mod( omega - tau_grid, omega );   
   
    Pext_vals = interp1(sP, HP(:,2), taus_to_eval, 'pchip');
    Pext_vals = max(0, min(1, Pext_vals));
    
    Pmat(ii, :) = 1 - Pext_vals;
    
end

% ----- contour plot -----

[TauMesh, KfMesh] = meshgrid(tau_grid,  alpha_values);
contourf(TauMesh, KfMesh, Pmat, 25, 'LineColor','none');
colorbar;
xlabel('t_0');
ylabel('${\alpha}$','Interpreter','latex');
colormap(jet);

