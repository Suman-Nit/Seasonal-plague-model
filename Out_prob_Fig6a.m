clc; close all;clear;


K_R   = 250000;       
N_H   = 10000;        
r_R   = 5/365;        
p     = 0.65;         
d_R   = 0.00055;      
m_R   = 0.0547;       
g_R   = 0.06;         
mu_R  = 0.03;         

a     = 3/250000;     
r_F   = 0.0084;       
d_H   = 0.00011;     
r_H   = 0.00011;      
m_H   = 0.1;          
g_H   = 0.34;         
omega_1 = 160;        

sigma1 = 0.150994; sigma2 = 0.033235; k_f=0.019684; alpha =  0.64195;  
beta_R_bar = 0.53057; beta_H_bar = 0.0097422;  K_F_bar = 3; d_F_bar = 0.20569; 
 
Tper = 365;

beta_R = @(t) beta_R_bar * (1 + sigma1 * cos(2*pi*t/Tper));
beta_H = @(t) beta_H_bar * (1 + sigma2 * cos(2*pi*t/Tper));
d_Ffun = @(t) d_F_bar    * (1 + alpha  * cos(2*pi*(t - omega_1)/Tper));
K_Ffun = @(t) K_F_bar * (1 + k_f).^(sin(2*pi*t/Tper));  


g1 = @(x1,x2,x3,tau) ((d_R+m_R) + (d_R+m_R*(1-g_R)).*K_Ffun(tau).*x1.*x2) ./ ...
                     ((d_R+m_R) + (d_R+m_R*(1-g_R)).*K_Ffun(tau));
g2 = @(x1,x2,x3,tau) ( beta_R(tau).*(1-exp(-a*K_R)).*x1.*x2 + d_Ffun(tau) + (1-exp(-a*K_R)) ...
                     + beta_H(tau).*N_H.*exp(-a*K_R).*x2.*x3 ) ...
                     ./ ( beta_R(tau).*(1-exp(-a*K_R)) + d_Ffun(tau) + (1-exp(-a*K_R)) ...
                     + beta_H(tau).*N_H.*exp(-a*K_R) );
g3 = @(x1,x2,x3,tau) 1;

dHdt = @(s,H,tabs) [ - ((d_R+m_R) + (d_R+m_R*(1-g_R))*K_Ffun(tabs - s)) * (H(1) - g1(H(1),H(2),H(3), tabs - s)); ...
                     - ( beta_R(tabs - s).*(1-exp(-a*K_R)) + d_Ffun(tabs - s) + (1-exp(-a*K_R)) ...
                       + beta_H(tabs - s).*N_H.*exp(-a*K_R) ) * (H(2) - g2(H(1),H(2),H(3), tabs - s)); ...
                     - (d_H + m_H) * (H(3) - g3(H(1),H(2),H(3), tabs - s)) ];

tabs = Tper;                
num_points = 10000;
s_span = linspace(0, 3*Tper, num_points);
H0 = [0; 0; 0];
opts = odeset('RelTol',1e-7,'AbsTol',1e-9);
[~, Hsol] = ode45(@(s,H) dHdt(s,H,tabs), s_span, H0, opts);

tol = 1e-9;
idxP = s_span >= 2*Tper - tol & s_span <= 3*Tper + tol;
sP   = s_span(idxP) - 2*Tper;       
HP   = Hsol(idxP,:);

H1_per = @(s) interp1(sP, HP(:,1), mod(s,Tper), 'pchip', 'extrap');
H2_per = @(s) interp1(sP, HP(:,2), mod(s,Tper), 'pchip', 'extrap');
H3_per = @(s) interp1(sP, HP(:,3), mod(s,Tper), 'pchip', 'extrap');

tau_grid = linspace(0, Tper, 1000);
Pext_1_0_0 = arrayfun(@(tau) H1_per(Tper - tau), tau_grid);   
Pext_0_1_0 = arrayfun(@(tau) H2_per(Tper - tau), tau_grid);   
Pext_0_0_1 = arrayfun(@(tau) H3_per(Tper - tau), tau_grid);   

Pout_1_0_0 = 1 - Pext_1_0_0;
Pout_0_1_0 = 1 - Pext_0_1_0;  
Pout_0_0_1 = 1 - Pext_0_0_1;

tau_values = linspace(0, 365, 10);   
tots       = 10000;                  
M_threshold = 150;                  
t_max      = 365;                   
rng(42);

Pout_MC = zeros(size(tau_values));

for k = 1:numel(tau_values)
    t0 = tau_values(k);
    n_outbreak = 0;

    for rep = 1:tots
        I_R = 1; R_R = 0; S_R = K_R - I_R - R_R;
        F   = 0;                      
        Nf  = K_Ffun(t0);            
        I_H = 0; R_H = 0; S_H = N_H - I_H - R_H;

        t = 0;                        
        sumI = I_R + F + I_H;

        while (sumI > 0) && (sumI < M_threshold) && (t < t_max)
            
            t_curr = t0 + t;

            T_R = S_R + I_R + R_R;
            if T_R <= 0
                break;
            end

            betaR = beta_R(t_curr);
            betaH = beta_H(t_curr);
            dF    = d_Ffun(t_curr);
            KF    = K_Ffun(t_curr);

            A = zeros(20,1);
            A(1)  = r_R*S_R*(1 - (T_R/K_R));                                    
            A(2)  = r_R*R_R*(1 - p);                                             
            A(3)  = d_R*S_R;                                                     
            A(4)  = betaR * (S_R/T_R) * F * (1 - exp(-a*T_R));                   
            A(5)  = F * (1 - exp(-a*T_R));                                       
            A(6)  = dF * F;                                                      
            A(7)  = (d_R + m_R*(1 - g_R)) * I_R;                                 
            A(8)  = r_R*R_R*(p - T_R/K_R);                                       
            A(9)  = m_R*g_R*I_R;                                                 
            A(10) = d_R*R_R;                                                     
            A(11) = r_F*Nf*(1 - Nf/KF);                                          
            A(12) = (F/T_R) * (1 - exp(-a*T_R));                                
            A(13) = (d_R + m_R*(1 - g_R)) * I_R * Nf;                            
            A(14) = r_H*(S_H + R_H);                                             
            A(15) = d_H*S_H;                                                     
            A(16) = betaH * N_H * F * exp(-a*T_R);                               
            A(17) = (d_H + m_H*(1 - g_H)) * I_H;                                 
            A(18) = m_H*g_H*I_H;                                                 
            A(19) = d_H*R_H;                                                    
            A(20) = m_H*(1 - g_H) * I_H;                                         

            B = sum(A);
            if B <= 0
                break;  
            end

            eps2 = rand;
            
              t  = t + 0.005;

            % Select event
            r = eps2 * B;
            c = 0; ev = 0;
            for ii = 1:20
                c = c + A(ii);
                if r <= c
                    ev = ii; break;
                end
            end
            if ev == 0, ev = 20; end
            
            switch ev
                case 1   
                    S_R = S_R + 1;
                case 2   
                    S_R = S_R + 1;
                case 3   
                    S_R = max(S_R - 1,0);
                case 4   
                    S_R = max(S_R - 1,0);  I_R = I_R + 1;
                case 5   
                    F   = max(F - 1,0);
                case 6   
                    F   = max(F - 1,0);
                case 7   
                    I_R = max(I_R - 1,0);
                case 8   
                    R_R = R_R + 1;
                case 9   
                    I_R = max(I_R - 1,0);  R_R = R_R + 1;
                case 10  
                    R_R = max(R_R - 1,0);
                case 11  
                    Nf  = Nf + 1;
                case 12  
                    Nf = Nf + 1;
                case 13  
                    F   = F + 1;
                case 14  
                    S_H = S_H + 1;
                case 15  
                    S_H = max(S_H - 1,0);
                case 16  
                    S_H = max(S_H - 1,0);  I_H = I_H + 1;
                case 17  
                    I_H = max(I_H - 1,0);
                case 18  
                    I_H = max(I_H - 1,0);  R_H = R_H + 1;
                case 19  
                    R_H = max(R_H - 1,0);
                case 20  
                    
            end

            sumI = I_R + F + I_H;
        end

        
        if sumI >= M_threshold
            n_outbreak = n_outbreak + 1;     
        end
    end

    Pout_MC(k) = n_outbreak / tots;
end

 ----------------
figure; hold on; box on;
plot(tau_grid, Pout_1_0_0, 'k-', 'LineWidth', 2.5, 'DisplayName','Analytical');
plot(tau_values, Pout_MC, 'bo', 'LineWidth', 2, 'MarkerSize', 6, ...
    'DisplayName','Numerical');

xlabel('Introduction time, t_0 (days)');
ylabel('Outbreak probability');
title('Outbreak probability vs introduction time (single infectious flea)');
legend('Location','best');
ylim([0 1]);

