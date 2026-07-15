% prior_plot.m
% Independent script to visualize the prior distribution
clear; clc; close all;

% 1. Define Prior Parameters
T_min = 0.1;
T_max = 2.0;
alpha = -5;

% 2. Create the prior model using the existing utility class
prior_model = prior_utils.create_prior_model(T_min, T_max, alpha);

% 3. Generate a temperature grid and evaluate the PDF
T_grid = linspace(0, T_max + 0.5, 1000); % Extended grid to show the boundaries
prior_pdf = prior_utils.evaluate_prior_pdf(T_grid, prior_model);

% 4. Plot the results
figure('Name', 'Prior Distribution', 'Position', [100, 100, 600, 400]);
plot(T_grid, prior_pdf, 'b-', 'LineWidth', 2);
hold on;
xline(T_min, 'r--', 'T_{min}', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');
xline(T_max, 'r--', 'T_{max}', 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');

% Formatting
title(sprintf('Prior Distribution (T_{min} = %.1f, T_{max} = %.1f, \\alpha = %.1f)', T_min, T_max, alpha));
xlabel('Temperature (T)');
ylabel('Probability Density p(T)');
grid on;
box on;
xlim([0, T_max + 0.5]);
set(gca, 'FontSize', 12);