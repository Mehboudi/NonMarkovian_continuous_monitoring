function [hypotheses, T_est, rel_err] = update_bayesian_filter(hypotheses, W_true, T_true, sys_vars, sys_matrices, params, dt, estimation_metric)
    % UPDATE_BAYESIAN_FILTER Progresses all hypotheses and returns Bayesian estimate
    
    % Get process noise for the vector of hypothesized temperatures
    flags.calc_Q = true;
    mats = evaluate_T_matrices(hypotheses.Tvec, sys_vars, sys_matrices, flags);
    Q_hypotheses = mats.Q;
    
    F = sys_matrices.F;
    
    % Predict Step (Vectorized across hypotheses)
    sigma_bar = pagemtimes(pagemtimes(F, hypotheses.sigma), F.') + Q_hypotheses;
    d_bar = F * hypotheses.d;
    
    mean_W = sqrt(params.lambda_meas * dt) * d_bar(3, :);
    var_W = 0.5 + params.lambda_meas * dt * (squeeze(sigma_bar(3, 3, :))' - 0.5);
    
    % Likelihood evaluation
    likelihoods = (1./sqrt(2*pi*var_W)) .* exp(-(W_true - mean_W).^2 ./ (2*var_W));
    
    % Extract kalman blocks
    C_mat = squeeze(sigma_bar(1:2, 3:4, :));       
    sigma_2 = squeeze(sigma_bar(3:4, 3:4, :)); 
    
    e1 = [1; 0];
    K_mode1 = squeeze(pagemtimes(C_mat, e1));             
    K_mode2 = squeeze(pagemtimes(sigma_2 - 0.5*eye(2), e1)); 
    K_unscaled = [K_mode1; K_mode2];                           
    K = K_unscaled ./ var_W;
    
    innovation = W_true - mean_W; 
    
    % Update Step
    hypotheses.d = d_bar + innovation * sqrt(params.lambda_meas * dt) .* K;
    
    K_reshaped = reshape(K, [4, 1, length(hypotheses.Tvec)]); 
    var_W_reshaped = reshape(var_W, [1, 1, length(hypotheses.Tvec)]);
    KVK = pagemtimes(K_reshaped .* var_W_reshaped, pagetranspose(K_reshaped)); 
    hypotheses.sigma = sigma_bar - params.lambda_meas * dt * KVK;
    
    hypotheses.weights = hypotheses.weights .* likelihoods;
    hypotheses.weights = prior_utils.normalize_probabilities(hypotheses.weights);

    % Estimator Calculation based on chosen metric
    if strcmp(estimation_metric, 'ERMSE')
        mean_inv_T = sum(hypotheses.weights .* (1./hypotheses.Tvec));
        mean_inv_T_sq = sum(hypotheses.weights .* (1./(hypotheses.Tvec.^2)));
        T_est = mean_inv_T / mean_inv_T_sq;
        rel_err = ((T_est - T_true) / T_true)^2;
    elseif strcmp(estimation_metric, 'EMSLE')
        T_est = exp(sum(hypotheses.weights .* log(hypotheses.Tvec)));
        rel_err = log(T_est / T_true)^2;
    end
end