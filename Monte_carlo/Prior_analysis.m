% prior_analysis.m
% Independent script to visualize the prior distribution and its information content
clear; clc; close all;

% Define Parameters
T_min = 0.1;
T_max = 2;
T_grid = linspace(T_min, T_max, 5000); % High resolution for accurate integration

figure('Name', 'Prior Distribution Analysis', 'Position', [100, 100, 1100, 450]);

%% 1. Plot the Prior Distribution (Example for alpha = -5)
alpha_example = 0;
prior_model_ex = prior_utils.create_prior_model(T_min, T_max, alpha_example);
prior_pdf_ex = prior_utils.evaluate_prior_pdf(T_grid, prior_model_ex);

subplot(1, 2, 1);
plot(T_grid, prior_pdf_ex, 'b-', 'LineWidth', 2);
hold on;
xline(T_min, 'r--', 'T_{min}', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');
xline(T_max, 'r--', 'T_{max}', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');

title(sprintf('Prior PDF (\\alpha = %.1f)', alpha_example));
xlabel('Temperature (T)');
ylabel('Probability Density p(T)');
grid on; box on;
xlim([0, T_max + 0.2]);
set(gca, 'FontSize', 12);

%% 2. Plot the Prior Information vs alpha
alpha_vec = linspace(-10, 10, 100);
prior_info_vec = zeros(size(alpha_vec));

for i = 1:length(alpha_vec)
    % Create model and evaluate PDF and its derivative
    pm = prior_utils.create_prior_model(T_min, T_max, alpha_vec(i));
    p = prior_utils.evaluate_prior_pdf(T_grid, pm);
    dp = prior_utils.evaluate_prior_pdf_derivative(T_grid, pm);
    
    % Calculate EMSLE prior information: \int [T*p'(T) + p(T)]^2 / p(T) dT
    % We use a threshold mask to prevent division by zero at the boundaries
    p_floor = 1e-12 * max(p);
    mask = p > p_floor;
    
    integrand = zeros(size(T_grid));
    numerator = (T_grid(mask) .* dp(mask) + p(mask)).^2;
    integrand(mask) = numerator ./ p(mask);
    
    prior_info_vec(i) = trapz(T_grid, integrand);
end

subplot(1, 2, 2);
plot(alpha_vec, prior_info_vec, 'k-', 'LineWidth', 2);
title('Prior Information vs. Shape Parameter (\alpha)');
xlabel('\alpha');
ylabel('Prior Information Q(T)');
grid on; box on;
set(gca, 'FontSize', 12);