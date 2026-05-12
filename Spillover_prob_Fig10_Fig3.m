clc; close all; clear;

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

% --- Seasonal parameters
sigma1 = 0.150994; 
sigma2 = 0.033235; 
k_f    = 0.019684; 
alpha  = 0.64195;  

beta_R_bar = 0.53057; 
beta_H_bar = 0.0097422;  
K_F_bar    = 3; 
d_F_bar    = 0.20569; 

Tper = 365;

% Seasonal functions
beta_R = @(t) beta_R_bar * (1 + sigma1 * cos(2*pi*t/Tper));
beta_H = @(t) beta_H_bar * (1 + sigma2 * cos(2*pi*t/Tper));
d_Ffun = @(t) d_F_bar    * (1 + alpha  * cos(2*pi*(t - omega_1)/Tper));
K_Ffun = @(t) K_F_bar * (1 + k_f).^(sin(2*pi*t/Tper));

%% ---------------- Analytical spillover probability ----------------
g1 = @(x1,x2,tau) ((d_R+m_R) + (d_R+m_R*(1-g_R)).*K_Ffun(tau).*x1.*x2) ./ ...
                  ((d_R+m_R) + (d_R+m_R*(1-g_R)).*K_Ffun(tau));

g2 = @(x1,x2,tau) ...
    ( beta_R(tau).*(1-exp(-a*K_R)).*x1.*x2 + ...
      d_Ffun(tau) + ...
      (1-exp(-a*K_R)) ) ./ ...
    ( beta_R(tau).*(1-exp(-a*K_R)) + ...
      d_Ffun(tau) + ...
      (1-exp(-a*K_R)) + ...
      beta_H(tau).*N_H.*exp(-a*K_R) );

dHdt = @(s,H,tabs) [ ...
    - ((d_R+m_R) + (d_R+m_R*(1-g_R))*K_Ffun(tabs-s)) .* ...
      (H(1) - g1(H(1),H(2),tabs-s)); ...
    - ( beta_R(tabs-s).*(1-exp(-a*K_R)) + ...
        d_Ffun(tabs-s) + ...
        (1-exp(-a*K_R)) + ...
        beta_H(tabs-s).*N_H.*exp(-a*K_R) ) .* ...
      (H(2) - g2(H(1),H(2),tabs-s)) ...
    ];

tabs = Tper;
num_points = 10000;
s_span = linspace(0,3*Tper,num_points);
H0 = [0;0];

opts = odeset('RelTol',1e-7,'AbsTol',1e-9);
[~,Hsol] = ode45(@(s,H) dHdt(s,H,tabs), s_span, H0, opts);

idxP = s_span >= 2*Tper & s_span <= 3*Tper;
sP   = s_span(idxP) - 2*Tper;
HP   = Hsol(idxP,:);

H1_per = @(s) interp1(sP,HP(:,1),mod(s,Tper),'pchip','extrap');
H2_per = @(s) interp1(sP,HP(:,2),mod(s,Tper),'pchip','extrap');

tau_grid = linspace(0,Tper,1000);

Pext_1_0 = arrayfun(@(tau) H1_per(Tper-tau), tau_grid);
Pext_0_1 = arrayfun(@(tau) H2_per(Tper-tau), tau_grid);

Pspill_analytic = 1 - Pext_0_1;

%% ---------------- CTMC Monte Carlo spillover simulation ----------------
tau_values   = linspace(0,365,10);
tots         = 10000;
t_max        = 365;

p_spill_sim = zeros(size(tau_values));

for k = 1:numel(tau_values)

    t0 = tau_values(k);
    spill_count = 0;

    for rep = 1:tots

        %% Initial state: one infectious flea
        I_R = 0;
        R_R = 0;
        S_R = K_R;

        F   = 1;
        Nf  = ceil(K_Ffun(t0));

        I_H = 0;
        R_H = 0;
        S_H = N_H;

        t = 0;

        while (F > 0) && (I_H == 0) && (t < t_max)

            T_R = S_R + I_R + R_R;

            if T_R <= 0
                break;
            end

            t_curr = t0 + t;

            betaR = beta_R(t_curr);
            betaH = beta_H(t_curr);
            dF    = d_Ffun(t_curr);
            KF    = K_Ffun(t_curr);

            %% Event rates
            A = zeros(16,1);

            A(1)  = r_R*S_R*(1 - T_R/K_R);
            A(2)  = r_R*R_R*(1-p);
            A(3)  = d_R*S_R;
            A(4)  = betaR*(S_R/T_R)*F*(1-exp(-a*T_R));
            A(5)  = F*(1-exp(-a*T_R));
            A(6)  = dF*F;
            A(7)  = (d_R+m_R*(1-g_R))*I_R;
            A(8)  = r_R*R_R*(p-T_R/K_R);
            A(9)  = m_R*g_R*I_R;
            A(10) = d_R*R_R;
            A(11) = r_F*Nf*(1-Nf/KF);
            A(12) = (F/T_R)*(1-exp(-a*T_R));
            A(13) = (d_R+m_R*(1-g_R))*I_R*Nf;
            A(14) = r_H*(S_H+R_H);
            A(15) = d_H*S_H;
            A(16) = betaH*S_H*F*exp(-a*T_R);

            B = sum(A);

            if B <= 0
                break;
            end

            %% Gillespie event selection
            r = rand * B;
            t = t + 0.005;

            c = 0;
            ev = 0;

            for ii = 1:16
                c = c + A(ii);
                if r <= c
                    ev = ii;
                    break;
                end
            end

            if ev == 0
                ev = 16;
            end

            %% State updates
            switch ev

                case 1
                    S_R = S_R + 1;

                case 2
                    S_R = S_R + 1;

                case 3
                    S_R = max(S_R - 1,0);

                case 4
                    S_R = max(S_R - 1,0);
                    I_R = I_R + 1;

                case 5
                    F = max(F - 1,0);

                case 6
                    F = max(F - 1,0);

                case 7
                    I_R = max(I_R - 1,0);

                case 8
                    R_R = R_R + 1;

                case 9
                    I_R = max(I_R - 1,0);
                    R_R = R_R + 1;

                case 10
                    R_R = max(R_R - 1,0);

                case 11
                    Nf = Nf + 1;

                case 12
                    F  = max(F - 1,0);
                    Nf = Nf + 1;

                case 13
                    F = F + 1;

                case 14
                    S_H = S_H + 1;

                case 15
                    S_H = max(S_H - 1,0);

                case 16
                    S_H = max(S_H - 1,0);
                    I_H = I_H + 1;

            end

        end

        %% Spillover classification
        if I_H > 0
            spill_count = spill_count + 1;
        end

    end

    p_spill_sim(k) = spill_count / tots;

end
beta_R = @(t) beta_R_bar * (1 + sigma1 * cos(2 * pi * t/365));
beta_H = @(t) beta_H_bar * (1 + sigma2 * cos(2 * pi * t/365));
d_F = @(t) d_F_bar * (1 + alpha * cos(2 * pi * (t - omega_1)/365));
K_F1 = @(t) K_F_bar * (1 + k_f).^(sin(2 * pi * t/365));
t = linspace(0, 365, 500); 
beta_R_values = beta_R(t); 
beta_H_values = beta_H(t);
d_F_values = d_F(t);
K_F_values = K_F1(t);

figure(1);
yyaxis left
plot(t, beta_R_values); ylabel('\beta_R(t)'); ylim([0.4,0.7])
yyaxis right
plot(t, beta_H_values); ylabel('\beta_H(t)'); ylim([0.006 0.011])
xlim([0 365]);
xlabel('t')

figure(2);
yyaxis left
plot(t, d_F_values); ylabel('d_F(t)'); ylim([0,0.5]);
yyaxis right
plot(t, K_F_values); ylabel('K_F(t)'); ylim([2.6 3.2]);
xlim([0 365]);
xlabel('t')

%% ---------------- Plot ----------------
figure(3);
hold on; box on;

plot(tau_grid, Pspill_analytic, 'k-', 'LineWidth', 2.5, ...
    'DisplayName','Analytical');

plot(tau_values, p_spill_sim, 'bo', ...
    'MarkerSize', 7, ...
    'LineWidth', 2, ...
    'DisplayName','Numerical');

xlabel('t_0 (day)');
ylabel('Spillover Probability');
title('Spillover probability vs introduction time');
legend('Location','best');

xlim([0 365]);
ylim([0 1]);