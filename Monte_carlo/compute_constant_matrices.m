function sys_matrices = compute_constant_matrices(params, sys_vars, dt)
    % COMPUTE_CONSTANT_MATRICES Builds A, F, D2, and Q basis tensors that 
    % are independent of temperature.
    
    w_p = sys_vars.w_p;
    w_m = sys_vars.w_m;
    theta = sys_vars.theta;
    J_p = sys_vars.J_p;
    J_m = sys_vars.J_m;
    lam_diss = sys_vars.lam_diss;
    
    % Constant drift coefficients
    cc1 = -sin(2*theta)/(4*w_p*w_m)*( J_m*w_p*cos(theta)^2 - J_p*w_m*sin(theta)^2 );
    dd1 = -sin(2*theta)^2/(8*w_p*w_m)*(J_m*w_p + J_p*w_m);
    dd2 = -1/(2*w_p*w_m)*(J_m*w_p*cos(theta)^4 + J_p*w_m*sin(theta)^4);
    
    A1 = [dd2, 1, cc1, 0;
          -sys_vars.w_rc^2, dd2, lam_diss, cc1;
          cc1, 0, dd1, 1;
          lam_diss, cc1, -sys_vars.w_s^2, dd1];
          
    A2 = zeros(4); A2(3:4, 3:4) = -params.lambda_meas / 2 * eye(2);
    
    D2 = zeros(4); D2(3:4, 3:4) = params.lambda_meas * 0.5 * eye(2);
    
    A = A1 + A2;
    F = expm(A * dt);
    
    % Process noise basis definition
    E_basis = zeros(4, 4, 6);
    E_basis(1,1,1) = 1; E_basis(2,2,2) = 1; 
    E_basis(3,3,3) = 1; E_basis(4,4,4) = 1;
    E_basis(1,3,5) = 1; E_basis(3,1,5) = 1;
    E_basis(2,4,6) = 1; E_basis(4,2,6) = 1;
    
    Q_basis = zeros(4, 4, 6);
    for ib = 1:6
        Q_basis(:,:,ib) = finite_horizon_process_noise(A, F, E_basis(:,:,ib));
    end
    Q_D2 = finite_horizon_process_noise(A, F, D2);

    % Store in struct
    sys_matrices = struct();
    sys_matrices.A = A;
    sys_matrices.F = F;
    sys_matrices.D2 = D2;
    sys_matrices.Q_basis = Q_basis;
    sys_matrices.Q_D2 = Q_D2;
    
    % Precompute measurement matrices for BCRB ARE equations
    sys_matrices.C_meas_full = [0, 0, 1, 0; 0, 0, 0, 1]; 
    sys_matrices.V_inv = [2, 0; 0, 0];
    sys_matrices.V_inv_sqrt = [sqrt(2), 0; 0, 0];
end

% --- Helper Function ---
function Q = finite_horizon_process_noise(A, F, D)
    C = D - F * D * F.';
    Q = lyap(A, A.', C);
    Q = 0.5 * (Q + Q.');
end