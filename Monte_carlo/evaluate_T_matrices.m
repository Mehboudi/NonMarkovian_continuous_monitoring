function out = evaluate_T_matrices(T, sys_vars, sys_matrices, flags)
    % EVALUATE_T_MATRICES Centralizes the temperature-dependent Math block.
    % Computes diffusion D(T), D'(T), and process noise Q(T) dynamically.
    
    w_p = sys_vars.w_p; w_m = sys_vars.w_m; theta = sys_vars.theta;
    J_p = sys_vars.J_p; J_m = sys_vars.J_m;
    
    N_p = sys_vars.N_bosonic(w_p, T);
    N_m = sys_vars.N_bosonic(w_m, T);
    
    % Compute shared coefficients
    c3 = 2*w_p^2*w_m^2;
    c1 = sin(2*theta)/4*(J_m*(2*N_m + 1)*cos(theta)^2 - J_p*(2*N_p + 1)*sin(theta)^2);
    c2 = sin(2*theta)/(2*c3)*(J_m*w_p^2*(2*N_m + 1)*cos(theta)^2 - J_p*w_m^2*(2*N_p + 1)*sin(theta)^2);
    d1 = sin(2*theta)^2/(4*c3)*(J_m*w_p^2*(2*N_m + 1) + J_p*w_m^2*(2*N_p + 1));
    d2 = sin(2*theta)^2/(8)*(J_m*(2*N_m + 1) + J_p*(2*N_p + 1));
    d3 = 1/c3 * (J_m*w_p^2*(2*N_m + 1)*cos(theta)^4 + J_p*w_m^2*(2*N_p + 1)*sin(theta)^4);
    d4 = 1/2 * (J_m*(2*N_m + 1)*cos(theta)^4 + J_p*(2*N_p + 1)*sin(theta)^4);
    
    out = struct();
    
    % 1. Diffusion Matrix (D)
    if isfield(flags, 'calc_D') && flags.calc_D
        out.D = sys_matrices.D2 + diag([d3, d4, d1, d2]);
        out.D(1,3) = c2; out.D(3,1) = c2;
        out.D(2,4) = c1; out.D(4,2) = c1;
    end
    
    % 2. Derivative of Diffusion Matrix (D')
    if isfield(flags, 'calc_D_prime') && flags.calc_D_prime
       dN_p = sys_vars.dN_dT(w_p, T);
       dN_m = sys_vars.dN_dT(w_m, T);
       
       dc1_dT = sin(2*theta)/4*(J_m*(2*dN_m)*cos(theta)^2 - J_p*(2*dN_p)*sin(theta)^2);
       dc2_dT = sin(2*theta)/(2*c3)*(J_m*w_p^2*(2*dN_m)*cos(theta)^2 - J_p*w_m^2*(2*dN_p)*sin(theta)^2);
       dd1_dT = sin(2*theta)^2/(4*c3)*(J_m*w_p^2*(2*dN_m) + J_p*w_m^2*(2*dN_p));
       dd2_dT = sin(2*theta)^2/(8)*(J_m*(2*dN_m) + J_p*(2*dN_p));
       dd3_dT = 1/c3 * (J_m*w_p^2*(2*dN_m)*cos(theta)^4 + J_p*w_m^2*(2*dN_p)*sin(theta)^4);
       dd4_dT = 1/2 * (J_m*(2*dN_m)*cos(theta)^4 + J_p*(2*dN_p)*sin(theta)^4);
       
       out.D_prime = diag([dd3_dT, dd4_dT, dd1_dT, dd2_dT]);
       out.D_prime(1,3) = dc2_dT; out.D_prime(3,1) = dc2_dT;
       out.D_prime(2,4) = dc1_dT; out.D_prime(4,2) = dc1_dT;
    end
    
    % 3. Process Noise Tensor (Q)
    if isfield(flags, 'calc_Q') && flags.calc_Q
        N_hyp = length(T);
        Q = repmat(sys_matrices.Q_D2, 1, 1, N_hyp);
        
        % NOTE: Maps directly to E_basis 1 to 6 layout created in compute_constant_matrices
        coeffs = [d3; d4; d1; d2; c2; c1]; 
        
        for ib = 1:6
            Q = Q + reshape(coeffs(ib, :), 1, 1, N_hyp) .* sys_matrices.Q_basis(:,:,ib);
        end
        
        % Return 2D matrix if scalar T, else 3D tensor
        out.Q = squeeze(Q); 
        if N_hyp == 1
             out.Q = reshape(out.Q, 4, 4);
        end
    end
end