function [final_res, best_KT, best_RT, success] = solve_shooting_2sector_inverted(p, Policy, Guess, Targets, Settings)
    if nargin < 5, Settings = struct(); end
    if ~isfield(Settings, 'tol'), Settings.tol = 1e-3; end
    if ~isfield(Settings, 'max_iter'), Settings.max_iter = 50; end
    if ~isfield(Settings, 'T_periods'), Settings.T_periods = 100; end
    if ~isfield(Settings, 'perturb_K'), Settings.perturb_K = 1e-10; end 
    if ~isfield(Settings, 'perturb_R'), Settings.perturb_R = 1e-10; end 
    if ~isfield(Settings, 'min_lambda'), Settings.min_lambda = 1e-8; end 
    
    T = Settings.T_periods;
    
    x_curr = [Guess.KT; Guess.RT];
    y_target = [Targets.K0; Targets.R0];
    
    fprintf('Executing 2-sector shooting algorithm...\n');
    
    success = false;
    final_res = [];
    best_KT = x_curr(1);
    best_RT = x_curr(2);
    
    fprintf('Verifying initial foothold...\n');
    [x_curr, valid_start] = ensure_valid_foothold_inverted(x_curr, p, T, Policy, 0.05); 
    
    if ~valid_start
        fprintf('Error: Invalid initial guess. Manually adjust Guess.KT.\n');
        return;
    end
    fprintf('Initial guess confirmed: KT=%.4f, RT=%.4f\n', x_curr(1), x_curr(2));
    
    for iter = 1:Settings.max_iter
        [y_curr, res_curr, valid_base, fail_msg] = run_simulation_wrapper_2sector(x_curr(1), x_curr(2), p, T, Policy);
        
        if ~valid_base
            fprintf('Iter %d: Baseline failure (%s). Attempting recovery...\n', iter, fail_msg);
            [x_curr, recovered] = ensure_valid_foothold_inverted(x_curr, p, T, Policy, 0.001); 
            if ~recovered
                fprintf('Recovery failed. Terminating.\n'); break;
            else
                 [y_curr, res_curr, ~, ~] = run_simulation_wrapper_2sector(x_curr(1), x_curr(2), p, T, Policy);
            end
        end
        
        error_vec = y_curr - y_target;
        err_norm = norm(error_vec);
        
        fprintf('Iter %2d: KT=%9.4f, RT=%9.4f -> Err=%.2e\n', iter, x_curr(1), x_curr(2), err_norm);
        
        if err_norm < Settings.tol
            fprintf('Convergence achieved.\n');
            success = true;
            final_res = res_curr;
            best_KT = x_curr(1);
            best_RT = x_curr(2);
            return;
        end
        
        J = zeros(2,2);
        
        dk = Settings.perturb_K * max(1, x_curr(1));
        [dk_col, found_k] = compute_gradient_safe(@(k) run_simulation_wrapper_2sector(k, x_curr(2), p, T, Policy), ...
                                                  x_curr(1), y_curr, dk);
        dr = Settings.perturb_R * max(1, x_curr(2));
        [dr_col, found_r] = compute_gradient_safe(@(r) run_simulation_wrapper_2sector(x_curr(1), r, p, T, Policy), ...
                                                  x_curr(2), y_curr, dr);
                                              
       if ~found_k || ~found_r
            fprintf('Gradient computation failed. Reducing perturbation step.\n');
            Settings.perturb_K = Settings.perturb_K * 0.5;
            Settings.perturb_R = Settings.perturb_R * 0.5;
            continue;
        end
        
        J(:,1) = dk_col; J(:,2) = dr_col;
        
        fprintf('    [Jacobian] Det = %11.4e | Rcond = %11.4e\n', det(J), rcond(J));
        
        if rcond(J) < 1e-12
            fprintf('Singular Jacobian detected. Using gradient descent.\n');
            delta_x = J' * error_vec * 0.01;
        else
            delta_x = J \ error_vec;
        end
        
        lambda = 1.0; 
        found_safe_step = false;
        
        while lambda > Settings.min_lambda 
            x_try = x_curr - lambda * delta_x;
            
            if x_try(1) < 10 || x_try(2) < 1 
                lambda = lambda * 0.5; continue;
            end
            
            [y_try, ~, valid_try, ~] = run_simulation_wrapper_2sector(x_try(1), x_try(2), p, T, Policy);
            
            if valid_try
                new_err = norm(y_try - y_target);
                
                allow_backdoor = (lambda < 0.1) && (err_norm > 0.5);
                
                if new_err < err_norm || allow_backdoor
                    x_curr = x_try;
                    found_safe_step = true;
                    fprintf('    -> Step accepted (lambda=%.1e): Error=%.2e\n', lambda, new_err);
                    break;
                else
                    lambda = lambda * 0.5;
                end
            else
                lambda = lambda * 0.5;
            end
        end
        
        if ~found_safe_step
            fprintf('Line search failed. Current error: %.2e\n', err_norm);
            
            if err_norm < 1e-3
                fprintf('Micro-scale error bound reached. Terminating successfully.\n');
                success = true;
                final_res = res_curr;
                best_KT = x_curr(1);
                best_RT = x_curr(2);
                return;
            end
            
            fprintf('Applying micro-perturbation...\n');
            x_curr = x_curr .* [1.0001; 1]; 
            Settings.perturb_K = Settings.perturb_K * 0.5;
        end
    end
    
    final_res = res_curr;
    best_KT = x_curr(1);
    best_RT = x_curr(2);
end

function [x_safe, success] = ensure_valid_foothold_inverted(x_start, p, T, Policy, max_adjust_ratio)
    x_curr = x_start;
    success = false;
    step_ratio = 1e-4;  
    
    for attempt = 1:50
        [~, ~, valid, msg] = run_simulation_wrapper_2sector(x_curr(1), x_curr(2), p, T, Policy);
        
        if valid
            success = true;
            if attempt > 1
                fprintf('Safe point found after %d adjustments: KT=%.7f\n', attempt-1, x_curr(1));
            end
            x_safe = x_curr;
            return;
        end
        
        if abs(x_curr(1) - x_start(1)) / x_start(1) > max_adjust_ratio
            fprintf('Adjustment exceeded %.1f%% limit. Aborting.\n', max_adjust_ratio*100);
            break;
        end
        
        if contains(msg, 'Depletion') || contains(msg, 'Capital') || contains(msg, 'Net Investment')
            fprintf('Attempt %d: Capital depletion -> Increasing KT (%.7f -> %.7f)\n', attempt, x_curr(1), x_curr(1)*(1+step_ratio));
            x_curr(1) = x_curr(1) * (1 + step_ratio);
        elseif contains(msg, 'Crowding') || contains(msg, 'Melt-down') || contains(msg, 'Investment')
            fprintf('Attempt %d: Investment crowding out -> Decreasing KT (%.7f -> %.7f)\n', attempt, x_curr(1), x_curr(1)*(1-step_ratio));
            x_curr(1) = x_curr(1) * (1 - step_ratio);
        else
            fprintf('Attempt %d: Unknown failure (%s) -> Applying random perturbation\n', attempt, msg);
            x_curr = x_curr .* (1 + step_ratio * randn(2,1));
        end
        
        step_ratio = step_ratio * 0.5; 
    end
    x_safe = x_start; 
end

function [grad, success] = compute_gradient_safe(func, x_val, y_base, delta)
    success = false;
    grad = [0; 0];
    
    for k = 1:10
        [y_plus, ~, valid_plus, ~] = func(x_val + delta);
        if valid_plus
            grad = (y_plus - y_base) / delta;
            success = true; return;
        end
        
        [y_minus, ~, valid_minus, ~] = func(x_val - delta);
        if valid_minus
            grad = (y_base - y_minus) / delta;
            success = true; return;
        end
        delta = delta * 0.1;
    end
end

function [y, res, is_valid, msg] = run_simulation_wrapper_2sector(KT, RT, p, T, Policy)
    msg = '';
    try
        [~, term] = evalc('solve_terminal_ABGP_2sector(KT, RT, p, T, Policy)');
    catch
        y = [nan; nan]; res = []; is_valid = false; msg = 'ABGP Crash'; return;
    end
    
    final_eq.r_next = term.r; final_eq.h_next = term.h;
    final_eq.K_next = term.K; final_eq.R_next = term.R;
    final_eq.E_next = term.E; final_eq.I_next = term.I; 
    
    try
        [~, res, debug_info] = evalc('solve_full_backward_path_2sector(final_eq, p, T, Policy)');
        if strcmp(debug_info.status, 'fail')
            y = [nan; nan]; is_valid = false;
            msg = debug_info.msg; 
            return;
        end
    catch
        y = [nan; nan]; res = []; is_valid = false; msg = 'Solver Exception'; return;
    end
    
    if any(isnan(res.K)) || res.K(1) <= 0 || res.R(1) <= 0
        y = [nan; nan]; is_valid = false; msg = 'Invalid End State';
    else
        y = [res.K(1); res.R(1)]; is_valid = true; msg = 'OK';
    end
end