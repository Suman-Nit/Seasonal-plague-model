
clc; clear; close all;

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
 
Tper = 365;             

beta_R = @(t) beta_R_bar * (1 + sigma1 * cos(2*pi*t/Tper));
beta_H = @(t) beta_H_bar * (1 + sigma2 * cos(2*pi*t/Tper));
d_Ffun = @(t) d_F_bar    * (1 + alpha  * cos(2*pi*(t - omega_1)/Tper));
K_Ffun = @(t) K_F_bar * (1 + k_f)^(sin(2*pi*t/Tper));


t0      = 0;            
T_hor   = 100;          
I0_rat  = 1;            
F0_flea = 1;            
I0_hum  = 0;            
N0_flea = 3;   


paths = cell(5,1);

for i = 1:5
    [t_abs, IR, Fvec, IH] = simulate_ctmc_path( ...
        t0, T_hor, I0_rat, F0_flea, I0_hum, N0_flea, ...
        K_R, N_H, r_R, p, d_R, m_R, g_R, a, r_F, d_H, r_H, m_H, g_H, ...
        beta_R, beta_H, d_Ffun, K_Ffun);
    paths{i} = struct('t', t_abs, 'IR', IR, 'F', Fvec, 'IH', IH);
end

tspan = linspace(0, T_hor, 1000);

K0 = [K_R - I0_rat, I0_rat, 0, N0_flea, F0_flea, N_H, 0, 0, 0];

dydt = @(t, K) [
    r_R * K(1) * (1 - (K(1)+K(2)+K(3))/K_R) + r_R * K(3) * (1 - p) ...
        - d_R * K(1) - beta_R(t) * K(1) * K(5) * (1 - exp(-a * (K(1)+K(2)+K(3)))) / max(K(1)+K(2)+K(3), 1);

    beta_R(t) * K(1) * K(5) * (1 - exp(-a * (K(1)+K(2)+K(3)))) / max(K(1)+K(2)+K(3), 1) ...
        - (d_R + m_R) * K(2);
    
    r_R * K(3) * (p - (K(1)+K(2)+K(3))/K_R) + m_R * g_R * K(2) - d_R * K(3);
    
    r_F * K(4) * (1 - K(4) / K_Ffun(t)) + K(5) * (1 - exp(-a * (K(1)+K(2)+K(3)))) / max(K(1)+K(2)+K(3), 1);
    
    (d_R + m_R * (1 - g_R)) * K(2) * K(4) - d_Ffun(t) * K(5) ...
        - K(5) * (1 - exp(-a * (K(1)+K(2)+K(3))));
    
    r_H * (K(6) + K(8)) - d_H * K(6) ...
        - beta_H(t) * K(6) * K(5) * exp(-a * (K(1)+K(2)+K(3)));
    
    beta_H(t) * K(6) * K(5) * exp(-a * (K(1)+K(2)+K(3))) - (d_H + m_H) * K(7);
    

    m_H * g_H * K(7) - d_H * K(8);
    
    m_H * (1 - g_H) * K(7)
];
opts = odeset('RelTol',1e-10,'AbsTol',1e-9);
[ts, K] = ode45(dydt, tspan, K0, opts);
I_R_ode = K(:,2);
F_ode = K(:,5);
I_H_ode = K(:,7);



figure(3);hold on; box on;
for i = 1:5
    plot(paths{i}.t, paths{i}.F, '-', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('CTMC path %d', i));
end
plot(ts, F_ode, 'k-', 'LineWidth', 2.5, 'DisplayName','ODE');
 xlabel('t(day)'); ylabel('Infectiousfleas,F');
% title('Infectious humans: CTMC paths vs Deterministic ODE');
% legend('Location','best'); xlim([0 T_hor]);




function [t_abs, IR, Fvec, IH] = simulate_ctmc_path( ...
    t0, T_hor, I0_rat, F0_flea, I0_hum, N0_flea, ...
    K_R, N_H, r_R, p, d_R, m_R, g_R, a, r_F, d_H, r_H, m_H, g_H, ...
    beta_R, beta_H, d_Ffun, K_Ffun)

    % Initial state
    I_R = I0_rat;  F = F0_flea;  I_H = I0_hum;
    R_R = 0;       S_R = K_R - I_R - R_R;
    R_H = 0;       S_H = N_H - I_H - R_H;
    Nf  = N0_flea;

    % Storage
    t = 0; t_max_path = T_hor;
    t_vec = 0;
    IR = I_R;
    Fvec = F;
    IH = I_H;

    while t < t_max_path
        t_curr = t0 + t;

        T_R = S_R + I_R + R_R;
        if T_R <= 0, break; end

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
        if B <= 0, break; end

         u2 = rand;
        
        t_next = t + 0.005;
        if t_next > t_max_path
            break;
        end

        rsel = u2 * B; csum = cumsum(A);
        ev = find(csum >= rsel, 1);
        switch ev
            case 1,   S_R = S_R + 1;
            case 2,   S_R = S_R + 1;
            case 3,   S_R = max(S_R - 1,0);
            case 4,   S_R = max(S_R - 1,0);  I_R = I_R + 1;
            case 5,   F   = max(F - 1,0);
            case 6,   F   = max(F - 1,0);
            case 7,   I_R = max(I_R - 1,0);
            case 8,   R_R = R_R + 1;
            case 9,   I_R = max(I_R - 1,0);  R_R = R_R + 1;
            case 10,  R_R = max(R_R - 1,0);
            case 11,  Nf  = Nf + 1;
            case 12,  Nf = Nf + 1;
            case 13,  F   = F + 1;
            case 14,  S_H = S_H + 1;
            case 15,  S_H = max(S_H - 1,0);
            case 16,  S_H = max(S_H - 1,0);  I_H = I_H + 1;
            case 17,  I_H = max(I_H - 1,0);
            case 18,  I_H = max(I_H - 1,0);  R_H = R_H + 1;
            case 19,  R_H = max(R_H - 1,0);
            case 20  
        end

        
        t = t_next;
        t_vec(end+1,1) = t; 
        IR(end+1,1) = I_R; 
        Fvec(end+1,1) = F;
        IH(end+1,1) = I_H;

        if (I_R + F + I_H) == 0
            break;
        end
    end

    t_abs = t0 + t_vec;
end
