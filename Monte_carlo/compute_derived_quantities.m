function sys_vars = compute_derived_quantities(params)
    % COMPUTE_DERIVED_QUANTITIES Calculates frequencies, spectral densities,
    % and the bosonic occupation functions.

    sys_vars = struct();
    sys_vars.w_s_bare = 1.0;
    
    % Derived constants
    sys_vars.lam_diss = sqrt(params.alpha1 * params.alpha2 * params.gamma);
    sys_vars.w_rc = sqrt(params.alpha2 * params.gamma);
    sys_vars.dw_s = sys_vars.lam_diss / sys_vars.w_rc;
    sys_vars.w_s = sqrt(sys_vars.w_s_bare^2 + sys_vars.dw_s^2);
    
    % Normal modes
    term_p_m = sqrt(4 * sys_vars.lam_diss^2 + (sys_vars.w_rc^2 - sys_vars.w_s^2)^2);
    sys_vars.w_p = sqrt(sys_vars.w_s^2 + sys_vars.w_rc^2 + term_p_m) / sqrt(2);
    sys_vars.w_m = sqrt(sys_vars.w_s^2 + sys_vars.w_rc^2 - term_p_m) / sqrt(2);
    
    % Mixing angle
    term1 = -sys_vars.w_rc^2 + sys_vars.w_s^2 + term_p_m; %%%%%Careful with this
    term2 = 2 * term_p_m;
    sys_vars.theta = acos(sqrt(term1 / term2));
    
    % Spectral density functions
    J_rc = @(w) params.gamma * w .* exp(-abs(w) / params.Lambda_cutoff);
    sys_vars.J_p = J_rc(sys_vars.w_p);
    sys_vars.J_m = J_rc(sys_vars.w_m);
    
    % Bosonic occupation and its derivative
    sys_vars.N_bosonic = @(w, T) 1 ./ expm1(w ./ T);
    sys_vars.dN_dT = @(w, T) (w ./ T.^2) .* exp(w ./ T) ./ (expm1(w ./ T).^2);
end