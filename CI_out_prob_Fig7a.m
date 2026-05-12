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

tau_values = linspace(0, 365, 10);   
M_threshold = 150;                  
t_max      = 365;                   

Pout_MC = zeros(size(tau_values));

tots       = 5000;     
num_batches = 300;     

Pout_mean  = zeros(size(tau_values));
Pout_lower = zeros(size(tau_values));
Pout_upper = zeros(size(tau_values));

for k = 1:numel(tau_values)

    t0 = tau_values(k);
    out_prob_samples = zeros(num_batches,1);

    for b = 1:num_batches
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
                if B <= 0, break; end

                dt = -log(rand)/B;
                t  = t + dt;

                r = rand * B;
                c = cumsum(A);
                ev = find(r <= c,1);

                switch ev
                    case 1,  S_R = S_R + 1;
                    case 2,  S_R = S_R + 1;
                    case 3,  S_R = max(S_R - 1,0);
                    case 4,  S_R = max(S_R - 1,0); I_R = I_R + 1;
                    case 5,  F   = max(F - 1,0);
                    case 6,  F   = max(F - 1,0);
                    case 7,  I_R = max(I_R - 1,0);
                    case 8,  R_R = R_R + 1;
                    case 9,  I_R = max(I_R - 1,0); R_R = R_R + 1;
                    case 10, R_R = max(R_R - 1,0);
                    case 11, Nf  = Nf + 1;
                    case 12, Nf = Nf + 1;
                    case 13, F   = F + 1;
                    case 14, S_H = S_H + 1;
                    case 15, S_H = max(S_H - 1,0);
                    case 16, S_H = max(S_H - 1,0); I_H = I_H + 1;
                    case 17, I_H = max(I_H - 1,0);
                    case 18, I_H = max(I_H - 1,0); R_H = R_H + 1;
                    case 19, R_H = max(R_H - 1,0);
                end

                sumI = I_R + F + I_H;
            end

            if sumI >= M_threshold
                n_outbreak = n_outbreak + 1;
            end
        end

        out_prob_samples(b) = n_outbreak / tots;
    end

    Pout_mean(k)  = mean(out_prob_samples);
    CI = prctile(out_prob_samples,[2.5 97.5]);
    Pout_lower(k) = CI(1);
    Pout_upper(k) = CI(2);
end

figure; hold on; box on;

fill([tau_values fliplr(tau_values)], ...
     [Pout_lower fliplr(Pout_upper)], ...
     [0.85 0.85 0.85], 'EdgeColor','none');

plot(tau_values, Pout_mean, 'k-', 'LineWidth',2);

plot(tau_values, Pout_mean, 'bo','MarkerFaceColor','b');

xlabel('Introduction time, t_0 (days)');
ylabel('Outbreak probability');
title('Outbreak probability with 95% credible interval');
legend('95% CI','Mean (CTMC)','Location','best');
ylim([0 1]);
