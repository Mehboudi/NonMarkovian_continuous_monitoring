function [true_sys, W_true] = simulate_true_system(true_sys, T_true, sys_vars, sys_matrices, params, dt)
    % SIMULATE_TRUE_SYSTEM Propagates true state variables and generates W_true
    
    % Get Process Noise Matrix dynamically based on T_true
    flags.calc_Q = true;
    mats = evaluate_T_matrices(T_true, sys_vars, sys_matrices, flags);
    Q_true = mats.Q;
    
    F = sys_matrices.F;
    
    % Advance true system state
    sigma_bar_true = F * true_sys.sigma * F.' + Q_true;
    d_bar_true = F * true_sys.d;
    
    mean_W_true = sqrt(params.lambda_meas * dt) * d_bar_true(3);
    var_W_true = 0.5 + params.lambda_meas * dt * (sigma_bar_true(3, 3) - 0.5);
    
    % Generate measurement innovation
    W_true = mean_W_true + sqrt(var_W_true) * randn();
    
    % Apply measurement update
    C_true = sigma_bar_true(1:2, 3:4);
    sigma_2_true = sigma_bar_true(3:4, 3:4);
    K_true = [C_true ; sigma_2_true - 0.5*eye(2)] * [1;0] / var_W_true;
    
    true_sys.sigma = sigma_bar_true - params.lambda_meas * dt * (K_true * var_W_true * K_true');
    true_sys.d = d_bar_true + (W_true - mean_W_true) * sqrt(params.lambda_meas * dt) * K_true;
end