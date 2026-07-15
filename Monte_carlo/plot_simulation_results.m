function plot_simulation_results(sim_data, bcrb, estimation_metric)
    % PLOT_SIMULATION_RESULTS Renders three core figures from simulation data.

    % Generate density matrix for Figure 1
    n_bins = 100;
    n_edges = linspace(sim_data.T_min, sim_data.T_max, n_bins + 1);
    n_centers = (n_edges(1:end-1) + n_edges(2:end)) / 2;
    posterior_density = zeros(n_bins, sim_data.save_counter);
    for i = 1:sim_data.save_counter
        for j = 1:n_bins
            bin_indices = (sim_data.T_history_sparse(:, i) >= n_edges(j)) & (sim_data.T_history_sparse(:, i) < n_edges(j+1));
            if any(bin_indices)
                posterior_density(j, i) = sum(sim_data.weights_history_sparse(bin_indices, i));
            end
        end
    end

    % Figure 1: Posterior Distribution
    figure('Name', 'Bayesian Estimation of Temperature [Run 1]', 'Position', [100, 100, 1000, 600]);
    contourf(sim_data.t_sparse, n_centers, posterior_density, 40, 'LineColor', 'none'); hold on;
    plot(sim_data.t, ones(size(sim_data.t)) * sim_data.sampled_true_temperatures(1), 'r-', 'LineWidth', 3, 'DisplayName', 'True T (Run 1)');
    plot(sim_data.t, sim_data.all_estimator_histories(1, :), 'w--', 'LineWidth', 2.5, 'DisplayName', 'Bayesian T_{est}');
    title('Posterior Distribution of Temperature T vs. Time');
    xlabel('Time'); ylabel('T Hypothesis');
    ylim([sim_data.T_min, sim_data.T_max]); xlim([0, sim_data.t_final]);
    grid on; box on; legend('show', 'Location', 'northwest');
    h = colorbar; ylabel(h, 'Posterior Probability Density');
    colormap(jet); set(gca, 'FontSize', 12);

    % Figure 2: Record & Errors
    figure('Name', 'Estimator Performance', 'Position', [1150, 100, 800, 700]);
    subplot(2, 1, 1);
    plot(sim_data.t, sim_data.W_true_history, 'k-', 'LineWidth', 1);
    title('Simulated Measurement Outcome W(t) [Run 1]');
    xlabel('Time'); ylabel('Outcome W'); grid on; xlim([0, sim_data.t_final]); set(gca, 'FontSize', 12);
    
    subplot(2, 1, 2);
    loglog(sim_data.t(2:end), sim_data.avg_relative_error_history(2:end), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Avg Bayesian Error'); hold on;
    loglog(sim_data.t(2:end), bcrb.BCRB_curve, 'r--', 'LineWidth', 2, 'DisplayName', 'Bayesian CRB');
    loglog(sim_data.t(2:end), bcrb.unbiased_crb_curve, 'g-.', 'LineWidth', 2, 'DisplayName', 'Unbiased-Estimator CRB');
    if strcmp(estimation_metric, 'ERMSE')
        title(sprintf('Average Relative Squared Error (%d MC Runs)', sim_data.num_m));
        ylabel('Mean ((T_{est} - T_{true}) / T_{true})^2');
    else
        title(sprintf('Average EMSLE (%d MC Runs)', sim_data.num_m));
        ylabel('Mean (ln(T_{est} / T_{true}))^2');
    end
    xlabel('Time'); grid on; legend('show', 'Location', 'southwest'); xlim([sim_data.dt, sim_data.t_final]); set(gca, 'FontSize', 12);

    % Figure 3: True System State Evolution
    figure('Name', 'True System Evolution', 'Position', [100, 800, 1850, 800]);
    subplot(3, 1, 1);
    plot(sim_data.t, sim_data.true_sigma_history(1, :), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Var(x_1)'); hold on;
    plot(sim_data.t, sim_data.true_sigma_history(2, :), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Var(p_1)');
    plot(sim_data.t, sim_data.true_sigma_history(3, :), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Var(x_2)');
    plot(sim_data.t, sim_data.true_sigma_history(4, :), 'm-', 'LineWidth', 1.5, 'DisplayName', 'Var(p_2)');
    title('Variances [Run 1]'); xlabel('Time'); grid on; xlim([0, sim_data.t_final]); legend('show', 'Location', 'best');

    subplot(3, 1, 2);
    plot(sim_data.t, sim_data.true_sigma_history(5, :), 'LineWidth', 1.5, 'DisplayName', 'Cov(x_1, p_1)'); hold on;
    plot(sim_data.t, sim_data.true_sigma_history(6, :), 'LineWidth', 1.5, 'DisplayName', 'Cov(x_1, x_2)');
    plot(sim_data.t, sim_data.true_sigma_history(7, :), 'LineWidth', 1.5, 'DisplayName', 'Cov(x_1, p_2)');
    plot(sim_data.t, sim_data.true_sigma_history(8, :), 'LineWidth', 1.5, 'DisplayName', 'Cov(p_1, x_2)');
    plot(sim_data.t, sim_data.true_sigma_history(9, :), 'LineWidth', 1.5, 'DisplayName', 'Cov(p_1, p_2)');
    plot(sim_data.t, sim_data.true_sigma_history(10, :), 'LineWidth', 1.5, 'DisplayName', 'Cov(x_2, p_2)');
    title('Covariances [Run 1]'); xlabel('Time'); grid on; xlim([0, sim_data.t_final]); legend('show', 'Location', 'best');

    subplot(3, 1, 3);
    plot(sim_data.t, sim_data.true_d_history(1, :), 'r-', 'LineWidth', 1.5, 'DisplayName', '<x_1>'); hold on;
    plot(sim_data.t, sim_data.true_d_history(2, :), 'b-', 'LineWidth', 1.5, 'DisplayName', '<p_1>');
    plot(sim_data.t, sim_data.true_d_history(3, :), 'g-', 'LineWidth', 1.5, 'DisplayName', '<x_2>');
    plot(sim_data.t, sim_data.true_d_history(4, :), 'm-', 'LineWidth', 1.5, 'DisplayName', '<p_2>');
    title('Mean Quadratures [Run 1]'); xlabel('Time'); grid on; xlim([0, sim_data.t_final]); legend('show', 'Location', 'best');
end