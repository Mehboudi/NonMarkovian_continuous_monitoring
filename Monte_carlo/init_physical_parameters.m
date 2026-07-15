function params = init_physical_parameters()
    % INIT_PHYSICAL_PARAMETERS Initializes fixed system parameters.
    % Returns a struct 'params' containing physical and configuration constants.

    params = struct();
    params.eps = 1e-1;
    params.alpha1 = params.eps * 1e1;
    params.alpha2 = 1e1;
    params.gamma = 1;            
    params.Lambda_cutoff = 1000;
    params.lambda_meas = 1;      % Measurement strength
end