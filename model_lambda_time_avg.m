function dominant_eigenvalue = model_lambda_time_avg(lambda, tspan, beta_R_bar, sigma1, beta_H_bar, sigma2,...
    d_F_bar, alpha, K_F_bar, k_f, d_R, m_R, a, K_R, g_R, N_H, d_H, m_H, omega_1)
    
    beta_R = @(t) beta_R_bar * (1 + sigma1 * cos(2 * pi * t / 365));
    beta_H = @(t) beta_H_bar * (1 + sigma2 * cos(2 * pi * t / 365));
    d_F = @(t) d_F_bar * (1 + alpha * cos(2 * pi * (t  - omega_1)/ 365));
    K_F1 = @(t) K_F_bar * (1 + k_f)^(sin(2 * pi * t / 365));

    
    A = @(t) [-(d_R + m_R), beta_R(t) * (1 - exp(-a * K_R)) / lambda, 0; ...
              (d_R + m_R * (1 - g_R)) * K_F1(t) / lambda, -d_F(t) - (1 - exp(-a * K_R)), 0; ...
              0, beta_H(t) * N_H * exp(-a * K_R) / lambda, -(d_H + m_H)];

    
    n = size(A(0), 1); 
    I = eye(n); 
    y0 = I(:); 

    odefunc = @(t, y) reshape(A(t) * reshape(y, n, n), [], 1);

    [~, Y] = ode45(odefunc, tspan, y0);
 
    FinalMatrix = reshape(Y(end, :), n, n);
    
    eigenvalues = eig(FinalMatrix);

    dominant_eigenvalue = max(abs(eigenvalues));
end
