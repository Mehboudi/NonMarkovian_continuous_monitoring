classdef prior_utils
    % PRIOR_UTILS Groups helper functions for prior PDF calculation and sampling.
    % Defined as a classdef with static methods to satisfy grouping in one file 
    % while allowing direct calling from external scripts.
    
    methods (Static)
        function prior_model = create_prior_model(theta_min, theta_max, alpha)
            prior_model.theta_min = theta_min;
            prior_model.theta_max = theta_max;
            prior_model.alpha = alpha;
        end
        
        function p = evaluate_prior_pdf(theta, prior_model)
            L = prior_model.theta_max - prior_model.theta_min;
            phase = pi * (theta - prior_model.theta_min) ./ L;
            sin2_term = sin(phase).^2;
            
            if abs(prior_model.alpha) < 1e-10
                p = (2 / L) .* sin2_term;
            else
                k_alpha = exp(prior_model.alpha / 2) * besseli(0, prior_model.alpha / 2) - 1;
                if abs(k_alpha) < 1e-12
                    p = (2 / L) .* sin2_term;
                else
                    p = (exp(prior_model.alpha .* sin2_term) - 1) ./ (k_alpha * L);
                end
            end
            p(theta < prior_model.theta_min | theta > prior_model.theta_max) = 0;
            p = max(p, 0);
        end
        
        function dp = evaluate_prior_pdf_derivative(theta, prior_model)
            L = prior_model.theta_max - prior_model.theta_min;
            phase = pi * (theta - prior_model.theta_min) ./ L;
            sin2_term = sin(phase).^2;
            dsin2_dtheta = (pi / L) .* sin(2 * phase);
            
            if abs(prior_model.alpha) < 1e-10
                dp = (2 / L) .* dsin2_dtheta;
            else
                k_alpha = exp(prior_model.alpha / 2) * besseli(0, prior_model.alpha / 2) - 1;
                if abs(k_alpha) < 1e-12
                    dp = (2 / L) .* dsin2_dtheta;
                else
                    dp = (exp(prior_model.alpha .* sin2_term) .* prior_model.alpha .* dsin2_dtheta) ./ (k_alpha * L);
                end
            end
            dp(theta < prior_model.theta_min | theta > prior_model.theta_max) = 0;
        end
        
        function samples = sample_from_prior(prior_model, output_size)
            grid_size = 5000;
            theta_grid = linspace(prior_model.theta_min, prior_model.theta_max, grid_size);
            prior_grid = prior_utils.evaluate_prior_pdf(theta_grid, prior_model);
            
            cdf_grid = cumsum(prior_grid);
            cdf_grid = cdf_grid ./ cdf_grid(end);
            
            [cdf_unique, unique_idx] = unique(cdf_grid, 'stable');
            theta_unique = theta_grid(unique_idx);
            
            if numel(cdf_unique) < 2
                samples = unifrnd(prior_model.theta_min, prior_model.theta_max, output_size);
                return;
            end
            
            u = rand(output_size);
            samples = interp1(cdf_unique, theta_unique, u, 'linear', 'extrap');
            samples = min(max(samples, prior_model.theta_min), prior_model.theta_max);
        end
        
        function prob = normalize_probabilities(prob)
            total = sum(prob);
            if total <= 0 || ~isfinite(total)
                prob = ones(size(prob)) / numel(prob);
            else
                prob = prob / total;
            end
        end
    end
end