clear; clc; close all;

%% 1. Simulation Controls
estimation_metric = 'EMSLE'; 
T_true_nominal    = 1;
t_final           = 5e4;
dt                = 0.01;
N_hypotheses      = 200;
num_m             = 321;

T_min = 0.1;
T_max = 2.0;
prior_alpha = 0;

t = 0:dt:t_final;
N_steps = length(t) - 1;

%% 2. Initialization & Precomputation
params = init_physical_parameters();
sys_vars = compute_derived_quantities(params);
sys_matrices = compute_constant_matrices(params, sys_vars, dt);

prior_model = prior_utils.create_prior_model(T_min, T_max, prior_alpha);

%% 3. BCRB and Fisher Information 
bcrb = compute_bcrb(params, sys_vars, sys_matrices, prior_model, t_final, dt, estimation_metric);

%% 4. Simulation Execution 
fprintf('Starting %d Monte Carlo simulations...\n', num_m);

all_estimator_histories = zeros(num_m, N_steps + 1);
avg_relative_error_history = zeros(1, N_steps + 1);
sampled_true_temperatures = zeros(num_m, 1);

% Sparse Plot tracking arrays
plot_save_interval = max(1, floor(N_steps / 500));
num_saved_steps = floor(N_steps / plot_save_interval) + 1;
sim_data.T_history_sparse = zeros(N_hypotheses, num_saved_steps);
sim_data.weights_history_sparse = zeros(N_hypotheses, num_saved_steps);
sim_data.t_sparse = zeros(1, num_saved_steps);

% === RUN 1: Tracking Execution ===
fprintf('Run 1 / %d\n', num_m);
sampled_true_temperatures(1) = T_true_nominal;

% Init Hypotheses
hypotheses.Tvec = linspace(T_min, T_max, N_hypotheses);
hypotheses.weights = prior_utils.normalize_probabilities(prior_utils.evaluate_prior_pdf(hypotheses.Tvec, prior_model));
hypotheses.d = zeros(4, N_hypotheses);
hypotheses.sigma = zeros(4, 4, N_hypotheses);
for i1 = 1:N_hypotheses
    hypotheses.sigma(:,:,i1) = 0.5 * diag([coth(sys_vars.w_rc./(2*hypotheses.Tvec(i1)))*ones(1,2),ones(1,2)]);
end

% Init True System
true_sys.d = zeros(4, 1);
true_sys.sigma = 0.5 * diag([coth(sys_vars.w_rc/(2*T_true_nominal))*ones(1,2),ones(1,2)]);

estimator_history = zeros(1, N_steps + 1);
relative_error_history = zeros(1, N_steps + 1);
sim_data.W_true_history = zeros(1, N_steps + 1);
sim_data.true_d_history = zeros(4, N_steps + 1);
sim_data.true_sigma_history = zeros(10, N_steps + 1);

sim_data.T_history_sparse(:, 1) = hypotheses.Tvec';
sim_data.weights_history_sparse(:, 1) = hypotheses.weights';
sim_data.t_sparse(1) = t(1);
sim_data.true_d_history(:, 1) = true_sys.d;
sim_data.true_sigma_history(:, 1) = [diag(true_sys.sigma); true_sys.sigma(1,2); true_sys.sigma(1,3); true_sys.sigma(1,4); true_sys.sigma(2,3); true_sys.sigma(2,4); true_sys.sigma(3,4)];
save_counter = 1;

progress_interval = round(N_steps / 10);

for i = 1:N_steps
    if mod(i, progress_interval) == 0
        fprintf('  Run 1 Progress: %d%%\n', round(i / N_steps * 100));
    end

    [true_sys, W_true] = simulate_true_system(true_sys, T_true_nominal, sys_vars, sys_matrices, params, dt);
    
    [hypotheses, T_est, rel_err] = update_bayesian_filter(hypotheses, W_true, T_true_nominal, sys_vars, sys_matrices, params, dt, estimation_metric);
    
    sim_data.W_true_history(i+1) = W_true;
    sim_data.true_d_history(:, i+1) = true_sys.d;
    sim_data.true_sigma_history(:, i+1) = [diag(true_sys.sigma); true_sys.sigma(1,2); true_sys.sigma(1,3); true_sys.sigma(1,4); true_sys.sigma(2,3); true_sys.sigma(2,4); true_sys.sigma(3,4)];
    
    estimator_history(i+1) = T_est;
    relative_error_history(i+1) = rel_err;
    
    if mod(i, plot_save_interval) == 0
        save_counter = save_counter + 1;
        sim_data.T_history_sparse(:, save_counter) = hypotheses.Tvec';
        sim_data.weights_history_sparse(:, save_counter) = hypotheses.weights';
        sim_data.t_sparse(save_counter) = t(i+1);
    end
end
sim_data.save_counter = save_counter;
all_estimator_histories(1, :) = estimator_history;
avg_relative_error_history = avg_relative_error_history + relative_error_history;

% === RUN 2 to num_m: Parallel Execution ===
if num_m > 1
    parfor m = 2:num_m
        fprintf('Started Run %d / %d\n', m, num_m);
        
        T_true_loc = prior_utils.sample_from_prior(prior_model, [1, 1]);
        sampled_true_temperatures(m) = T_true_loc;
        
        true_sys_loc.d = zeros(4, 1);
        true_sys_loc.sigma = 0.5 * diag([coth(sys_vars.w_rc/(2*T_true_loc))*ones(1,2),ones(1,2)]);
        
        hyp_loc.Tvec = linspace(T_min, T_max, N_hypotheses);
        hyp_loc.weights = prior_utils.normalize_probabilities(prior_utils.evaluate_prior_pdf(hyp_loc.Tvec, prior_model));
        hyp_loc.d = zeros(4, N_hypotheses);
        hyp_loc.sigma = zeros(4, 4, N_hypotheses);
        for i1 = 1:N_hypotheses
            hyp_loc.sigma(:,:,i1) = 0.5 * diag([coth(sys_vars.w_rc./(2*hyp_loc.Tvec(i1)))*ones(1,2),ones(1,2)]);
        end
        
        estimator_hist_loc = zeros(1, N_steps + 1);
        rel_err_hist_loc = zeros(1, N_steps + 1);
        
        for i = 1:N_steps
            [true_sys_loc, W_true_loc] = simulate_true_system(true_sys_loc, T_true_loc, sys_vars, sys_matrices, params, dt);
            
            [hyp_loc, T_est, rel_err] = update_bayesian_filter(hyp_loc, W_true_loc, T_true_loc, sys_vars, sys_matrices, params, dt, estimation_metric);
            
            estimator_hist_loc(i+1) = T_est;
            rel_err_hist_loc(i+1) = rel_err;
        end
        all_estimator_histories(m, :) = estimator_hist_loc;
        avg_relative_error_history = avg_relative_error_history + rel_err_hist_loc;
    end
end

avg_relative_error_history = avg_relative_error_history / num_m;

%% 5. Plotting Wrap-Up
% Compile payload for the visualization script
sim_data.t = t;
sim_data.dt = dt;
sim_data.t_final = t_final;
sim_data.T_min = T_min;
sim_data.T_max = T_max;
sim_data.num_m = num_m;
sim_data.sampled_true_temperatures = sampled_true_temperatures;
sim_data.all_estimator_histories = all_estimator_histories;
sim_data.avg_relative_error_history = avg_relative_error_history;

plot_simulation_results(sim_data, bcrb, estimation_metric);
save('matlab.mat', '-v7.3')

% Save a lightweight payload dedicated to plotting/debugging.
plot_payload.W_true_history = sim_data.W_true_history;
plot_payload.T_history_sparse = sim_data.T_history_sparse;
plot_payload.weights_history_sparse = sim_data.weights_history_sparse;
plot_payload.t_sparse = sim_data.t_sparse;
plot_payload.t = t;
plot_payload.avg_relative_error_history = avg_relative_error_history;
plot_payload.all_estimator_history_run1 = all_estimator_histories(1, :);
plot_payload.sampled_true_temperature_run1 = sampled_true_temperatures(1);
plot_payload.T_min = T_min;
plot_payload.T_max = T_max;
plot_payload.dt = dt;
plot_payload.t_final = t_final;
plot_payload.num_m = num_m;
plot_payload.bcrb = bcrb;
save('matlab_plot_payload.mat', '-struct', 'plot_payload', '-v7.3');
