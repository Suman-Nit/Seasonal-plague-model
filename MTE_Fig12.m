
clear; close all; clc;

K_R   = 250000;         
N_H   = 10000;          
r_R   = 5/365;          
p     = 0.65;           
d_R   = 0.00055;        
m_R   = 0.0547;         
g_R   = 0.06;           

a     = 3/250000;       
r_F   = 0.0084;         
d_H   = 0.00011;        
r_H   = 0.00011;        
m_H   = 0.1;            
g_H   = 0.34;           
omega_1 = 160;          

sigma1 = 0.150994; sigma2 = 0.033235; k_f=0.019684; alpha =  0.64195; 
beta_R_bar = 0.53057; beta_H_bar = 0.0097422;  K_F_bar = 3; d_F_bar = 0.20569;


beta_R = @(t) beta_R_bar.*(1 + sigma1.*cos(2*pi.*t./365));
beta_H = @(t) beta_H_bar.*(1 + sigma2.*cos(2*pi.*t./365));
d_F    = @(t) d_F_bar.*(1 + alpha.*cos(2*pi.*(t - omega_1)./365));
K_F1   = @(t) K_F_bar .* (1 + k_f).^(sin(2*pi.*t./365));   


omega = 365;                
T_max = 1000;                
dt = 2.0;                   
t_grid = 0:dt:T_max;
opts_ivp = odeset('RelTol',1e-6,'AbsTol',1e-8);


g1 = @(x1, x2, x3, tau) ( (d_R + m_R) + (d_R + m_R*(1 - g_R)) .* K_F1(tau) .* x1 .* x2 ) ./ ...
    ( (d_R + m_R) + (d_R + m_R*(1 - g_R)) .* K_F1(tau) );

g2 = @(x1, x2, x3, tau) ( beta_R(tau).*(1 - exp(-a .* K_R)).*x1.*x2 + d_F(tau) + (1 - exp(-a .* K_R)) + beta_H(tau).*N_H.*exp(-a .* K_R).*x2.*x3 ) ./ ...
    ( beta_R(tau).*(1 - exp(-a .* K_R)) + d_F(tau) + (1 - exp(-a .* K_R)) + beta_H(tau).*N_H.*exp(-a .* K_R) );

g3 = @(x1, x2, x3, tau) 1;

dHds = @(s, H, t_param) [ ...
    - ( (d_R + m_R) + (d_R + m_R.*(1 - g_R)) .* K_F1(t_param - s) ) .* ( H(1) - g1(H(1), H(2), H(3), t_param - s) ); ...
    - ( ( beta_R(t_param - s).*(1 - exp(-a .* K_R)) + d_F(t_param - s) + (1 - exp(-a .* K_R)) + beta_H(t_param - s).*N_H.*exp(-a .* K_R) ) ) .* ...
    ( H(2) - g2(H(1), H(2), H(3), t_param - s) ); ...
    - ( (d_H + m_H) ) .* ( H(3) - g3(H(1), H(2), H(3), t_param - s) ) ...
    ];

%% 
t_fixed = omega;
num_points = 8000;
s_span = linspace(0, 3*omega, num_points);
H0 = [0; 0; 0];
opts_periodic = odeset('RelTol',1e-8,'AbsTol',1e-9);

[~, H_traj] = ode45(@(s,H) dHds(s, H, t_fixed), s_span, H0, opts_periodic);

idx_periodic = s_span >= 2*omega - 1e-9 & s_span <= 3*omega + 1e-9;
s_periodic = s_span(idx_periodic) - 2*omega;      
H_periodic = H_traj(idx_periodic, :);

H1_periodic = @(s) interp1(s_periodic, H_periodic(:,1), mod(s,omega), 'pchip', 'extrap');
H2_periodic = @(s) interp1(s_periodic, H_periodic(:,2), mod(s,omega), 'pchip', 'extrap');
H3_periodic = @(s) interp1(s_periodic, H_periodic(:,3), mod(s,omega), 'pchip', 'extrap');


n_points = 12;
tau_points = (1:n_points) .* (omega / n_points);

MTE_100 = zeros(size(tau_points)); SD_100 = zeros(size(tau_points));
MTE_010 = zeros(size(tau_points)); SD_010 = zeros(size(tau_points));
MTE_001 = zeros(size(tau_points)); SD_001 = zeros(size(tau_points));

fprintf('Starting exact MTE computation (modified model) for %d tau points, T_max=%.0f, dt=%.3g\n', ...
    length(tau_points), T_max, dt);

for k = 1:length(tau_points)
    tau = tau_points(k);
    Pext_100 = H1_periodic(omega - tau);
    Pext_010 = H2_periodic(omega - tau);
    Pext_001 = H3_periodic(omega - tau);
    
    P100_t = zeros(size(t_grid));
    P010_t = zeros(size(t_grid));
    P001_t = zeros(size(t_grid));
    
    for m = 1:length(t_grid)
        tt = t_grid(m);
        if tt == 0
            P100_t(m) = 0;
            P010_t(m) = 0;
            P001_t(m) = 0;
        else
            t_total = tau + tt;
            [~, Hs] = ode45(@(s,H) dHds(s,H,t_total), [0 tt], [0;0;0], opts_ivp);
            Hend = Hs(end, :);
            P100_t(m) = Hend(1);
            P010_t(m) = Hend(2);
            P001_t(m) = Hend(3);
        end
    end
    
    eps_d = 1e-12;
    denom100 = max(Pext_100, eps_d);
    denom010 = max(Pext_010, eps_d);
    denom001 = max(Pext_001, eps_d);
    
    integrand_100 = 1 - (P100_t ./ denom100);
    integrand_010 = 1 - (P010_t ./ denom010);
    integrand_001 = 1 - (P001_t ./ denom010);
    
    integrand_100(integrand_100 < 0) = 0;
    integrand_010(integrand_010 < 0) = 0;
    integrand_001(integrand_001 < 0) = 0;
    
    MTE_100(k) = trapz(t_grid, integrand_100);
    MTE_010(k) = trapz(t_grid, integrand_010);
    MTE_001(k) = trapz(t_grid, integrand_001);
    
    I1_100 = trapz(t_grid, 2 .* t_grid .* integrand_100);
    I1_010 = trapz(t_grid, 2 .* t_grid .* integrand_010);
    I1_001 = trapz(t_grid, 2 .* t_grid .* integrand_001);
    
    Var100 = max(I1_100 - MTE_100(k)^2, 0);
    Var010 = max(I1_010 - MTE_010(k)^2, 0);
    Var001 = max(I1_001 - MTE_001(k)^2, 0);
    
    SD_100(k) = sqrt(Var100);
    SD_010(k) = sqrt(Var010);
    SD_001(k) = sqrt(Var001);
    
    fprintf('tau=%.1f  Pext10=%.4f  MTE100=%.2f  MTE010=%.2f  MTE001=%.2f\n ', tau, Pext_100, MTE_100(k), MTE_010(k), MTE_001(k));
end

figure(1);
errorbar(tau_points, MTE_100, SD_100, 'o-','LineWidth',1.4,'MarkerFaceColor','k');
xlabel('t_0 (day)'); ylabel('Time (days)');

figure(2);
errorbar(tau_points, MTE_010, SD_010, 'o-','LineWidth',1.4,'MarkerFaceColor','k');
xlabel('t_0 (day)'); ylabel('Time (days)');

figure(3);
errorbar(tau_points, MTE_001, SD_001, 'o-','LineWidth',1.4,'MarkerFaceColor','k');
xlabel('t_0 (day)'); ylabel('Time(days)');


