function bcrb = compute_bcrb(params, sys_vars, sys_matrices, prior_model, t_final, dt, estimation_metric)
    % COMPUTE_BCRB Evaluates the Bayesian Cramer-Rao bound over the prior grid.

    fprintf('Precomputing Fisher Information and Bayesian CRB...\n');
    T_grid = linspace(prior_model.theta_min, prior_model.theta_max, 1000);
    F_rate_grid = zeros(size(T_grid));
    
    A = sys_matrices.A;
    lambda_meas = params.lambda_meas;
    V_inv = sys_matrices.V_inv;
    V_inv_sqrt = sys_matrices.V_inv_sqrt;
    C_meas_full = sys_matrices.C_meas_full;
    
    % Request flags for D and D' matrices
    flags.calc_D = true;
    flags.calc_D_prime = true;
    
    for idx = 1:length(T_grid)
        T_val = T_grid(idx);
        
        % Fetch modular matrices
        mats = evaluate_T_matrices(T_val, sys_vars, sys_matrices, flags);
        D_val = mats.D;
        D_prime = mats.D_prime;
        
        % --- ROBUST STEADY-STATE SOLVER (Replaces unstable Euler loop) ---
        % Construct matrices for the Algebraic Riccati Equation (ARE)
        W_mat = lambda_meas * (C_meas_full' * V_inv * C_meas_full);
        
        Sigma_env_full = zeros(4, 4);
        Sigma_env_full(3:4, 3:4) = 0.5 * eye(2);
        
        A_care = A + Sigma_env_full * W_mat;
        Q_care = D_val - Sigma_env_full * W_mat * Sigma_env_full;
        B_care = sqrtm(W_mat); % B*B' = W_mat
        
        % MATLAB's care(A, B, Q) solves: A'*X + X*A - X*B*B'*X + Q = 0
        % We need: A_care*X + X*A_care' - X*W_mat*X + Q_care = 0 
        % Passing A_care' makes the equations match exactly.
        sigma_ss = care(A_care', B_care, Q_care);
        % -----------------------------------------------------------------
        
        % Get steady state full K matrix
        C_ss = sigma_ss(1:2, 3:4);
        sigma_S_ss = sigma_ss(3:4, 3:4);
        K_ss_full = [C_ss; sigma_S_ss - 0.5*eye(2)];
        
        % A_tilde is mathematically guaranteed to be Hurwitz (stable) by care()
        A_tilde = A - lambda_meas * K_ss_full * V_inv * C_meas_full;
        sigma_prime = lyap(A_tilde, D_prime);
        
        C_prime = sigma_prime(1:2, 3:4);
        sigma_S_prime = sigma_prime(3:4, 3:4);
        K_prime_full = [C_prime; sigma_S_prime]; 
        
        A_prime = zeros(4, 4); 
        cal_A = [A, zeros(4,4); A_prime, A_tilde];
        L_V_inv_L = C_meas_full' * V_inv * C_meas_full;
        cal_Q = [zeros(4,4), zeros(4,4); zeros(4,4), L_V_inv_L];
        
        Z0 = [K_ss_full; K_prime_full] * V_inv_sqrt;
        cal_P = lyap(cal_A', cal_Q); 
        
        F_rate_grid(idx) = lambda_meas^2 * trace(Z0' * cal_P * Z0);
    end
    
    % Integration of Fisher Information and Prior Rate
    prior_density_grid = prior_utils.evaluate_prior_pdf(T_grid, prior_model);
    fi_integrand = prior_density_grid .* (T_grid.^2 .* F_rate_grid);
    fi_information = trapz(T_grid, fi_integrand);
    
    prior_density_derivative_grid = prior_utils.evaluate_prior_pdf_derivative(T_grid, prior_model);
    prior_integrand = zeros(size(T_grid));
    p_floor = 1e-12 * max(prior_density_grid);
    stable_mask = prior_density_grid > p_floor;
    
    metric_shift = strcmp(estimation_metric, 'EMSLE');
    numerator = (T_grid .* prior_density_derivative_grid + metric_shift * prior_density_grid).^2;
    prior_integrand(stable_mask) = numerator(stable_mask) ./ prior_density_grid(stable_mask);
    prior_information = trapz(T_grid, prior_integrand);
    
    % Save Curves
    t_vec = 0:dt:t_final;
    bcrb.BCRB_curve = 1 ./ (fi_information .* t_vec(2:end) + prior_information);
    
    unbiased_crb_constant = trapz(T_grid, prior_density_grid ./ max(T_grid.^2 .* F_rate_grid, realmin));
    bcrb.unbiased_crb_curve = unbiased_crb_constant ./ t_vec(2:end);
end