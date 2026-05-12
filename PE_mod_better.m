
niter = 80000;
burnin = 30000;
 thin = 5;
 nchains = 4;

% ---------------- Observed weekly deaths (given) ------------------------
obs_rat = [23,17,27,45,27,80,54,71,110,110,124,153,184,154,135,173,156,196,243,221,...
           218,184,166,173,168,107,39,26,37,33,10,5,12,5,6,9,7,5,4,7,19,18,13,13,4,8,8,3,4,6,3,5]';
obs_hum = [0,0,3,1,3,9,3,6,14,18,18,31,41,44,36,69,86,83,91,106,132,107,99,46,33,22,26,...
           17,12,6,9,10,6,2,1,4,2,1,2,2,1,2,0,0,1,2,0,0,2,0,0,0]';

n_obs = numel(obs_rat);
if numel(obs_hum) ~= n_obs, error('rat and human observation vectors must have same length'); end
t_obs = (7:7:7*(n_obs))';

% ---------------- fixed parameters & baseline initial state -------------
fixed.r_r = 5/365;      fixed.phi_d_F = 160;   fixed.d_H = 0.00011;
fixed.a   = 3/1035;      fixed.K_r     = 250000;fixed.p   = 0.65;
fixed.d_r = 0.00055;    fixed.m_r     = 0.0547;fixed.g_r = 0.06;
fixed.r_f = 0.0084;     fixed.g_h     = 0.34;  fixed.r_h = 0.00011;
fixed.m_h = 0.1;        fixed.T_period= 365;

y0_base = zeros(10,1);
y0_base(1)=1035; y0_base(2)=5; y0_base(3)=0; y0_base(4)=24; % rats
y0_base(5)=7;   y0_base(6)=0;                              % fleas
y0_base(7)=1500; y0_base(8)=8; y0_base(9)=0; y0_base(10)=0; % humans


% ---------------- prior bounds (user provided) -------------------------
prior.lb = [ 0.001, 0,   0,     0,   1,   0,   0,   0 ]';
prior.ub = [ 1,     1,   0.2,   1,   11,  1,   1,   1 ]';
lb = prior.lb(:)'; ub = prior.ub(:)';
d = numel(lb);


manual_init_4 = [
    0.01, 0.05, 0.01, 0.02, 3, 0.05, 0.05, 0.05;
    0.3,  0.2,  0.05, 0.1,  8, 0.3,  0.3,  0.3;
    0.8,  0.7,  0.15, 0.5, 10, 0.8,  0.7,  0.6;
    0.1,  0.9,  0.1,  0.3,  5, 0.2,  0.6,  0.4
];

manual_init_4 = max(repmat(lb,4,1), min(repmat(ub,4,1), manual_init_4));

eps_adapt = 1e-6;
adapt_interval = 100;      
t0_adapt = 5000;           
base_c_cov = 2.38^2 / d;
tuning_factor = 0.2;    
c_cov = base_c_cov * tuning_factor;

chains_all = cell(nchains,1);
logpost_all = cell(nchains,1);

acceptRates = zeros(nchains,1);


if isempty(gcp('nocreate')), parpool('local'); end

parfor ch = 1:nchains

    
    theta0 = manual_init_4(ch,:);
    y0_chain = y0_base;
   
    chain = zeros(niter, d);
    chain(1,:) = theta0;
    logpost_curr = logposterior(theta0, y0_chain, lb, ub, obs_rat, obs_hum, t_obs, fixed);
   
    logpost_chain = -Inf(niter,1);
    if isfinite(logpost_curr)
      logpost_chain(1) = logpost_curr;
    end

    
    param_sd = 0.1 * (ub - lb);
    C = diag(param_sd.^2);
    C_prop = c_cov * C;
    C_prop2 = C_prop * 0.2;

    accept_count = 0;
    accept_count_stage2 = 0;

    for it = 2:niter
        theta_curr = chain(it-1,:);
       
         z1 = mvnrnd(zeros(1, d), C_prop) + theta_curr;
         z1 = reflect_bounds(z1,lb,ub);
         
        logpost_z1 = logposterior(z1, y0_chain, lb, ub, obs_rat, obs_hum, t_obs, fixed);
       
        logalpha1 = logpost_z1 - logpost_curr;
        if isfinite(logpost_z1) && log(rand) < logalpha1
            chain(it,:) = z1;
            logpost_chain(it) = logpost_z1;
            logpost_curr = logpost_z1;
            accept_count = accept_count + 1;
        else
            
              
             z2 = mvnrnd(zeros(1, d), C_prop2) + theta_curr;
             z2 = reflect_bounds(z2,lb,ub);
            logpost_z2 = logposterior(z2, y0_chain, lb, ub, obs_rat, obs_hum, t_obs, fixed);

            if isfinite(logpost_z2)
                logalpha2 = logpost_z2 - logpost_curr;
                if log(rand) < logalpha2
                    chain(it,:) = z2;
                    logpost_chain(it) = logpost_z2;
                    logpost_curr = logpost_z2;
                    accept_count_stage2 = accept_count_stage2 + 1;
                else
                    chain(it,:) = theta_curr;
                    logpost_chain(it) = logpost_curr;
                end
            else
                chain(it,:) = theta_curr;
                logpost_chain(it) = logpost_curr;
            end
        end


            if it > t0_adapt && mod(it, adapt_interval) == 0
                
                i0 = max(2, floor(it*0.5));
                sample_mat = chain(i0:it, :);
                S = cov(sample_mat) + eps_adapt * eye(d);
                C_prop  = c_cov * S;
                C_prop2 = 0.2 * C_prop;   
            end

    end

    chains_all{ch} = chain;
    logpost_all{ch} = logpost_chain;
    acceptRates(ch) = (accept_count + accept_count_stage2) / (niter-1);
   
    fprintf('Chain %d finished. Acceptance Rate = %.4f\n', ch, acceptRates(ch));
end

posterior_samples_all = cell(nchains, 1);
for ch = 1:nchains
    idx = (burnin+1):thin:niter;
    posterior_samples_all{ch} = chains_all{ch}(idx, :);
end

nsamps_per_chain = cellfun(@(x) size(x,1), posterior_samples_all);
if isempty(nsamps_per_chain) || min(nsamps_per_chain) < 2
    warning('Not enough samples were generated after burn-in and thinning to calculate R-hat.');
    Rhat = NaN(1,d);
else
    nmin = min(nsamps_per_chain);
    if any(nsamps_per_chain ~= nmin)
        warning('Chains have different sample lengths; truncating to %d samples.', nmin);
        for ch = 1:nchains, posterior_samples_all{ch} = posterior_samples_all{ch}(1:nmin, :); end
    end
    n = nmin; m = nchains;
    S_tensor = zeros(n, d, m);
    for ch = 1:m, S_tensor(:,:,ch) = posterior_samples_all{ch}; end

    Rhat = zeros(1,d);
    for j = 1:d
        chain_means = squeeze(mean(S_tensor(:,j,:), 1));
        B = n * var(chain_means);
        W = mean(squeeze(var(S_tensor(:,j,:), 0, 1)));
        var_hat = ((n-1)/n)*W + (1/n)*B;
        Rhat(j) = sqrt(var_hat / W);
    end
end

param_names = {'beta_rr_base','A_beta_rr','beta_hh_base','A_beta_hh',...
               'K_f_base','A_K_f','d_F_base','A_d_F'};
fprintf('\n--- Gelman-Rubin R-hat per parameter ---\n');
for j = 1:d
    fprintf(' %12s: R-hat = %.4f\n', param_names{j}, Rhat(j));
end

all_post = vertcat(posterior_samples_all{:});
if ~isempty(all_post)
    theta_hat = mean(all_post,1);
else
    theta_hat = NaN(1,d);
end
result.Rhat = Rhat;
result.theta_hat = theta_hat;
result.all_post = all_post;

if ~isempty(all_post) && all(isfinite(theta_hat))
    
    mu_post = mean(all_post, 1);
    lo = prctile(all_post, 2.5, 1);
    hi = prctile(all_post, 97.5, 1);

    fprintf('\n--- Posterior means and 95%% credible intervals ---\n');
    for j = 1:d
        fprintf(' %12s: mean = %.5g, 95%%%% CrI [%.5g, %.5g]\n', ...
            param_names{j}, mu_post(j), lo(j), hi(j));
    end

    
    [fit_rat_mu, fit_hum_mu] = simulate_weekly_deaths(t_obs, theta_hat, y0_base, fixed);

    
    figure(2);
    plot(t_obs, obs_rat);     % observed rats (weekly)
    hold on;
    plot(t_obs, fit_rat_mu, 'LineWidth', 2);             % fitted mean curve
    hold off; grid on; xlim([min(t_obs) max(t_obs)]);
    xlabel('Day'); ylabel('Weekly rat deaths');
    legend('Observed', 'Fitted (post. mean)', 'Location', 'best');
    title('Ratsmod');


     figure(3);
    plot(t_obs, obs_hum);     % observed humans (weekly)
    hold on;
    plot(t_obs, fit_hum_mu, 'LineWidth', 2);             % fitted mean curve
    hold off; grid on; xlim([min(t_obs) max(t_obs)]);
    xlabel('Day'); ylabel('Weekly human deaths');
    legend('Observed', 'Fitted (post. mean)', 'Location', 'best');
    title('Humansmod');
else
    warning('No posterior samples available or theta\_hat invalid; skipping fitted-curve and CI output.');
end

% % Plotting
% figure('Name', 'Trace Plotsmod');
% 
% for j=1:d
%     subplot(ceil(d/2),2,j);
% 
%      hold on;
%     for ch=1:m, plot(chains_all{ch}(:,j)); 
%     end
%     hold off;
%     title(param_names{j}, 'Interpreter', 'none');
% end
% 


function x = reflect_bounds(x, lb, ub)
    for i = 1:length(x)
        if x(i) < lb(i)
            x(i) = lb(i) + (lb(i) - x(i));
        elseif x(i) > ub(i)
            x(i) = ub(i) - (x(i) - ub(i));
        end
        x(i) = max(lb(i), min(ub(i), x(i)));
    end
end

function lp = logposterior(theta, y0_chain, lb, ub, obs_rat, obs_hum, t_obs, fixed)
    theta = theta(:)';
    if any(theta < lb) || any(theta > ub) || any(~isfinite(theta))
        lp = -Inf;
        return;
    end
    try
        [weekly_rat_mu, weekly_hum_mu] = simulate_weekly_deaths(t_obs, theta, y0_chain, fixed);
        ll_rat = sum(obs_rat .* log(weekly_rat_mu) - weekly_rat_mu);
        ll_hum = sum(obs_hum .* log(weekly_hum_mu) - weekly_hum_mu);
        lp = ll_rat + ll_hum;
        if ~isreal(lp), lp = -Inf; end
    catch
        lp = -Inf;
    end
end

function [weekly_rat, weekly_hum] = simulate_weekly_deaths(tt, theta, y0_in, fixed)
    t_max = max(tt);
    tspan = [0, t_max];
    opts = odeset('RelTol',1e-5,'AbsTol',1e-7, 'NonNegative', 1:10);
    [Tsol, Ysol] = ode45(@(t,y) plague_model_rhs(t,y,theta,fixed), tspan, y0_in, opts);
   
    D_r_at = interp1(Tsol, Ysol(:,4), tt);
    D_h_at = interp1(Tsol, Ysol(:,10), tt);
   
    weekly_rat = max(0, diff([0; D_r_at]));
    weekly_hum = max(0, diff([0; D_h_at]));
end

function dydt = plague_model_rhs(t,y,theta,fixed)
    S_r=y(1); I_r=y(2); R_r=y(3); D_r=y(4);
    N_f=y(5); F_f=y(6);
    S_h=y(7); I_h=y(8); R_h=y(9); D_h=y(10);
    T_r = S_r + I_r + R_r;

    beta_rr_base=theta(1); A_beta_rr=theta(2); beta_hh_base=theta(3); A_beta_hh=theta(4);
    K_f_base=theta(5); A_K_f=theta(6); d_F_base=theta(7); A_d_F=theta(8);

    r_r=fixed.r_r; phi_d_F=fixed.phi_d_F; d_H=fixed.d_H; a=fixed.a; K_r=fixed.K_r;
    p=fixed.p; d_r=fixed.d_r; m_r=fixed.m_r; g_r=fixed.g_r; r_f=fixed.r_f;
    g_h=fixed.g_h; r_h=fixed.r_h; m_h=fixed.m_h; T_period=fixed.T_period;

    Tr_safe = max(T_r, 1e-9);
    seasonal_term = cos(2*pi*t/T_period);
    beta_rr_t = beta_rr_base * (1 + A_beta_rr*seasonal_term);
    beta_hh_t = beta_hh_base * (1 + A_beta_hh*seasonal_term);
    d_F_t     = d_F_base     * (1 + A_d_F*cos(2*pi*(t - phi_d_F)/T_period));
    K_f_t     = K_f_base     * ((1 + A_K_f)^(sin(2*pi*t/T_period)));

    dS_r = r_r*S_r*(1 - T_r/K_r) + r_r*R_r*(1 - p) - d_r*S_r - beta_rr_t*S_r*F_f*(1 - exp(-a*T_r))/Tr_safe;
    dI_r = beta_rr_t*S_r*F_f*(1 - exp(-a*T_r))/Tr_safe - (d_r + m_r)*I_r;
    dR_r = r_r*R_r*(p - T_r/K_r) + m_r*g_r*I_r - d_r*R_r;
    dD_r = (1 - g_r)*m_r*I_r;
    dN_f = r_f*N_f*(1 - N_f/K_f_t) + F_f*(1 - exp(-a*T_r))/Tr_safe;
    dF_f = (d_r + m_r*(1 - g_r))*N_f*I_r - d_F_t*F_f - F_f*(1 - exp(-a*T_r));
    dS_h = r_h*(S_h + R_h) - d_H*S_h - beta_hh_t*S_h*F_f*(exp(-a*T_r));
    dI_h = beta_hh_t*S_h*F_f*(exp(-a*T_r)) - (d_H + m_h)*I_h;
    dR_h = m_h*g_h*I_h - d_H*R_h;
    dD_h = m_h*(1 - g_h)*I_h;

    dydt = [dS_r; dI_r; dR_r; dD_r; dN_f; dF_f; dS_h; dI_h; dR_h; dD_h];
end