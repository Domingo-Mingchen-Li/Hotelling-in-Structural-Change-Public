function [path_results, debug_info] = solve_full_backward_path_2sector(final_eq_vectors, params, T_periods, Policy)
    debug_info.status = 'ok';
    debug_info.fail_t = -1;
    debug_info.msg = 'Success';
    
    if nargin < 4
        Policy.tau_c_m_path = zeros(1, T_periods+1);
        Policy.tau_c_s_path = zeros(1, T_periods+1);
        Policy.tau_x_path = zeros(1, T_periods+1);
        Policy.tau_int_path = zeros(1, T_periods+1);
        Policy.tau_e_path = zeros(1, T_periods+1);
    end
    
    fprintf('Initiating 2-sector backward shooting algorithm...\n');
    
    num_points = T_periods + 1;
    fields = {'r','w','h','p_int','p_int_user',...
              'p_m','p_s',...
              'E','I','K','R','R_input',...
              'c_m','c_s',...
              'l_m','l_s','l_x','l_e',...
              'Y_GDP',...
              'A_m','A_s','A_x','A_e',...
              'P_star'}; 
    for i = 1:length(fields), res.(fields{i}) = nan(1, num_points); end
    
    idx_T = num_points; 
    
    res.K(idx_T) = final_eq_vectors.K_next;
    res.R(idx_T) = final_eq_vectors.R_next;
    res.r(idx_T) = final_eq_vectors.r_next;
    res.h(idx_T) = final_eq_vectors.h_next;
    
    if isfield(final_eq_vectors, 'E_next') && isfield(final_eq_vectors, 'I_next')
        res.E(idx_T) = final_eq_vectors.E_next;
        res.I(idx_T) = final_eq_vectors.I_next;
    else
        error('CRITICAL: E_next and I_next must be provided by ABGP solver.');
    end
    
    res.A_m(idx_T) = params.A0.m * (params.g_A_m ^ T_periods);
    res.A_s(idx_T) = params.A0.s * (params.g_A_s ^ T_periods);
    res.A_x(idx_T) = params.A0.x * (params.g_A_x ^ T_periods);
    res.A_e(idx_T) = params.A0.e * (params.g_A_e ^ T_periods);
    
    tech_T.A_m = res.A_m(idx_T); tech_T.A_s = res.A_s(idx_T); 
    tech_T.A_x = res.A_x(idx_T); tech_T.A_e = res.A_e(idx_T);
    
    Pol_T.tau_c_m = Policy.tau_c_m_path(idx_T);
    Pol_T.tau_c_s = Policy.tau_c_s_path(idx_T);
    Pol_T.tau_x   = Policy.tau_x_path(idx_T);
    Pol_T.tau_int = Policy.tau_int_path(idx_T);
    Pol_T.tau_e   = Policy.tau_e_path(idx_T);
    
    [w_T, f_prod_T, f_user_T, p_m_T, p_s_T] = local_reconstruct_prices(res.r(idx_T), res.h(idx_T), params, tech_T, Pol_T);
    
    res.w(idx_T) = w_T; 
    res.p_int(idx_T) = f_prod_T; res.p_int_user(idx_T) = f_user_T;
    res.p_m(idx_T) = p_m_T; res.p_s(idx_T) = p_s_T;
    
    p_m_g = p_m_T * (1 + Pol_T.tau_c_m);
    p_s_g = p_s_T * (1 + Pol_T.tau_c_s);
    res.P_star(idx_T) = (p_m_g^params.omega_m) * (p_s_g^params.omega_s);
    
    res.Y_GDP(idx_T) = res.E(idx_T) + res.I(idx_T) * (1 + Pol_T.tau_x);
    g_rates = get_growth_rates_2sector(params);
    res.R_input(idx_T) = res.R(idx_T) * (1 - g_rates.g_e);
    
    fprintf('Terminal state (T): r=%.4f, K=%.2f, A_m=%.2f\n', res.r(idx_T), res.K(idx_T), res.A_m(idx_T));
    
    current_guess_r = res.r(idx_T);
    
    for t_idx = idx_T : -1 : 2
        t_current = t_idx - 1; 
        t_future = t_idx;      
        
        future_state.r_next = res.r(t_future);
        future_state.h_next = res.h(t_future);
        future_state.E_next = res.E(t_future);     
        future_state.K_next = res.K(t_future);
        future_state.R_next = res.R(t_future);
        future_state.P_star_next = res.P_star(t_future);
        
        res.A_m(t_current) = res.A_m(t_future) / params.g_A_m;
        res.A_s(t_current) = res.A_s(t_future) / params.g_A_s;
        res.A_x(t_current) = res.A_x(t_future) / params.g_A_x;
        res.A_e(t_current) = res.A_e(t_future) / params.g_A_e;
        
        tech_levels.A_m = res.A_m(t_current);
        tech_levels.A_s = res.A_s(t_current);
        tech_levels.A_x = res.A_x(t_current);
        tech_levels.A_e = res.A_e(t_current);
        
        CurrentPol.tau_c_m = Policy.tau_c_m_path(t_current);
        CurrentPol.tau_c_s = Policy.tau_c_s_path(t_current);
        CurrentPol.tau_x   = Policy.tau_x_path(t_current);
        CurrentPol.tau_int = Policy.tau_int_path(t_current);
        CurrentPol.tau_e   = Policy.tau_e_path(t_current);
        CurrentPol.tau_x_next = Policy.tau_x_path(t_future); 
        
        guess.r = current_guess_r;
        [sol, success] = solve_backward_step_2sector(future_state, params, tech_levels, guess, CurrentPol, t_current);
        
        if ~success
            fprintf('Solver failed at t=%d. Running diagnostics...\n', t_current-1);
            
            [reason, details] = diagnose_failure(current_guess_r, future_state, params, tech_levels, CurrentPol);
            
            fprintf('Diagnostics: %s\n', reason);
            fprintf('Details: %s\n', details);
            
            debug_info.status = 'fail';
            debug_info.fail_t = t_current - 1;
            debug_info.msg = reason;
            break; 
        end
        
        threshold_I = params.delta * sol.K;
        if sol.I < threshold_I - 1e-4
            fprintf('Investment bound violated at t=%d (I=%.4f < %.4f).\n', t_current-1, sol.I, threshold_I);
            fprintf('Diagnostics: Capital decumulation detected.\n');
            
            res.K(1:t_current) = nan; 
            
            debug_info.status = 'fail';
            debug_info.fail_t = t_current - 1;
            debug_info.msg = 'Investment Meltdown (Decumulation)';
            break; 
        end
        
        res.r(t_current) = sol.r; res.w(t_current) = sol.w; res.h(t_current) = sol.h;
        res.p_int(t_current) = sol.p_int; res.p_int_user(t_current) = sol.p_int_user;
        res.p_m(t_current) = sol.p_m; res.p_s(t_current) = sol.p_s;
        
        res.E(t_current) = sol.E; 
        res.I(t_current) = sol.I;
        res.K(t_current) = sol.K; 
        res.R(t_current) = sol.R; 
        res.R_input(t_current) = sol.R_input;
        
        res.c_m(t_current) = sol.c_m; res.c_s(t_current) = sol.c_s;
        res.l_m(t_current) = sol.l_m; res.l_s(t_current) = sol.l_s;
        res.l_x(t_current) = sol.l_x; res.l_e(t_current) = sol.l_e;
        res.P_star(t_current) = sol.P_star;
        
        res.Y_GDP(t_current) = sol.E + sol.I * (1 + CurrentPol.tau_x);
        
        current_guess_r = sol.r;
    end
    
    path_results = res;
end

function [reason, detail_str] = diagnose_failure(r_guess, future, params, tech, Pol)
    h_net = future.h_next / (1 - params.delta + future.r_next);
    [w, f_prod, f_user, p_m, p_s] = local_reconstruct_prices(r_guess, h_net, params, tech, Pol);
    
    p_m_g = p_m * (1 + Pol.tau_c_m);
    p_s_g = p_s * (1 + Pol.tau_c_s);
    P_star = (p_m_g ^ params.omega_m) * (p_s_g ^ params.omega_s);
    
    term_future = (future.E_next^(params.epsilon - 1)) / (future.P_star_next^params.epsilon * (1 + Pol.tau_x_next));
    euler_rhs = params.beta * (1 - params.delta + future.r_next) * term_future;
    E = (euler_rhs * (P_star^params.epsilon) * (1 + Pol.tau_x)) ^ (1 / (params.epsilon - 1));
    
    exp_m = params.epsilon * params.omega_m + params.zeta_m;
    exp_s = params.epsilon * params.omega_s + params.zeta_s;
    Pi_star = (p_m_g ^ exp_m) * (p_s_g ^ exp_s);
    term_pigl = params.xi * (E^(1 - params.epsilon)) * Pi_star * (P_star^params.epsilon);
    
    c_m = (params.omega_m * E + params.zeta_m * term_pigl) / p_m_g;
    c_s = (params.omega_s * E + params.zeta_s * term_pigl) / p_s_g;
    
    val_L_m = (1 - params.alpha_m - params.beta_m) * p_m * c_m;
    val_L_s = (1 - params.alpha_s - params.beta_s) * p_s * c_s;
    L_cons = (val_L_m + val_L_s) / w;
    
    val_M_m = params.beta_m * p_m * c_m;
    val_M_s = params.beta_s * p_s * c_s;
    M_cons = (val_M_m + val_M_s) / f_user;
    
    gamma_x = 1 - params.alpha_x - params.beta_x;
    gamma_e = 1 - params.alpha_e - params.beta_e;
    
    term_A = gamma_e * f_prod / w;
    term_B = params.beta_x * w / (gamma_x * f_user);
    Const_LE = term_A * M_cons;
    Coeff_LE_LX = term_A * term_B;
    
    l_x_implied = (1 - L_cons - Const_LE) / (1 + Coeff_LE_LX);
    
    if l_x_implied < 0
        reason = 'Crowding Out';
        total_L_cons_chain = L_cons + Const_LE; 
        detail_str = sprintf('Total labor demand exceeds 1.0 (L_total=%.2f, L_cons=%.2f, L_e_cons=%.2f).', ...
            total_L_cons_chain, L_cons, Const_LE);
        return;
    end
    
    I_real = w * l_x_implied / gamma_x;
    K_supply_needed = (future.K_next - I_real) / (1 - params.delta);
    
    if K_supply_needed <= 0
        reason = 'Capital Depletion';
        detail_str = sprintf('Implied K(t)=%.2f < 0 for K(t+1)=%.2f (I=%.2f).', ...
            K_supply_needed, future.K_next, I_real);
        return;
    end
    
    reason = 'Numerical Instability';
    detail_str = 'L_cons and K_supply within bounds. Potential local minimum or complex roots in fzero.';
end

function [w, f_prod, f_user, p_m, p_s] = local_reconstruct_prices(r, h_net, p, tech, Pol)
    h_gross = h_net * (1 + Pol.tau_e);
    cost_h_log = log(h_gross);
    get_log_C = @(A, al, be, ga) -log(A) - al*log(al) - be*log(be) - ga*log(ga);
    
    C_x = get_log_C(tech.A_x, p.alpha_x, p.beta_x, 1-p.alpha_x-p.beta_x);
    C_e = get_log_C(tech.A_e, p.alpha_e, p.beta_e, 1-p.alpha_e-p.beta_e);
    
    gamma_x = 1 - p.alpha_x - p.beta_x;
    gamma_e = 1 - p.alpha_e - p.beta_e;
    wedge_int_log = log(1 + Pol.tau_int);
    
    term_const = C_x + p.beta_x * (C_e + wedge_int_log + p.beta_e * cost_h_log);
    term_r     = (p.alpha_x + p.beta_x * p.alpha_e) * log(r);
    coeff_w    = gamma_x + p.beta_x * gamma_e;
    
    w = exp(-(term_const + term_r) / coeff_w);
    
    f_prod = exp(C_e + p.alpha_e*log(r) + gamma_e*log(w) + p.beta_e*cost_h_log);
    f_user = f_prod * (1 + Pol.tau_int);
    
    calc_p = @(C, al, be, ga) exp(C + al*log(r) + ga*log(w) + be*log(f_user));
    C_m = get_log_C(tech.A_m, p.alpha_m, p.beta_m, 1-p.alpha_m-p.beta_m);
    p_m = calc_p(C_m, p.alpha_m, p.beta_m, 1-p.alpha_m-p.beta_m);
    C_s = get_log_C(tech.A_s, p.alpha_s, p.beta_s, 1-p.alpha_s-p.beta_s);
    p_s = calc_p(C_s, p.alpha_s, p.beta_s, 1-p.alpha_s-p.beta_s);
end

function visualize_2sector_results(res, T, Policy, params)
    if isnan(res.K(1))
        start_idx = find(~isnan(res.K), 1);
        if isempty(start_idx), return; end
    else
        start_idx = 1;
    end
    
    valid_idx = start_idx:(T+1);
    time = (start_idx-1):T;
    
    if length(valid_idx) < 2
        fprintf('Insufficient valid data points for visualization.\n');
        return;
    end
    
    figure('Name', '2-Sector Model Dynamics', 'Color', 'w', 'Position', [100, 100, 1400, 900]);
    
    subplot(3,3,1);
    yyaxis left; plot(time, res.K(valid_idx), 'b-', 'LineWidth', 2); ylabel('Capital K');
    yyaxis right; plot(time, res.R(valid_idx), 'g--', 'LineWidth', 2); ylabel('Resource R');
    title('1. State Variables'); grid on; xlim([0 T]);
    
    subplot(3,3,2);
    yyaxis left; plot(time, res.r(valid_idx), 'r-', 'LineWidth', 2); ylabel('Interest Rate r');
    yyaxis right; plot(time, res.h(valid_idx), 'k--', 'LineWidth', 2); ylabel('Net Resource Price h');
    title('2. Factor Prices'); grid on; xlim([0 T]);
    
    subplot(3,3,3);
    Val_m = res.p_m(valid_idx) .* (1 + Policy.tau_c_m_path(valid_idx)) .* res.c_m(valid_idx);
    Val_s = res.p_s(valid_idx) .* (1 + Policy.tau_c_s_path(valid_idx)) .* res.c_s(valid_idx);
    Total_C = Val_m + Val_s;
    
    plot(time, Val_m ./ Total_C, 'b-', 'LineWidth', 2); hold on;
    plot(time, Val_s ./ Total_C, 'r-', 'LineWidth', 2);
    legend('Manufacturing', 'Services');
    title('3. Expenditure Shares'); grid on; xlim([0 T]);
    
    subplot(3,3,4);
    area(time, [res.l_m(valid_idx)', res.l_s(valid_idx)', res.l_x(valid_idx)', res.l_e(valid_idx)']);
    legend('Manuf', 'Serv', 'Invest', 'Energy', 'Location', 'best');
    title('4. Labor Allocation'); xlim([0 T]); ylim([0 1]);
    
    subplot(3,3,5);
    inv_rate = res.I(valid_idx) ./ res.Y_GDP(valid_idx);
    plot(time, inv_rate, 'k-', 'LineWidth', 2);
    yline(0, 'r--');
    title('5. Investment Rate (I/Y)'); grid on; xlim([0 T]);
    
    subplot(3,3,6);
    plot(time, res.R_input(valid_idx), 'm-', 'LineWidth', 2);
    title('6. Extraction Flow'); grid on; xlim([0 T]);
    
    subplot(3,3,7);
    plot(time, res.p_m(valid_idx) ./ res.p_s(valid_idx), 'b-', 'LineWidth', 1.5);
    title('7. Relative Price (Pm / Ps)'); grid on; xlim([0 T]);
    
    subplot(3,3,8);
    plot(time, Policy.tau_e_path(valid_idx), 'r-', 'LineWidth', 2);
    title('8. Carbon Tax (\tau_e)'); grid on; xlim([0 T]);
    
    subplot(3,3,9);
    plot(time, Policy.tau_int_path(valid_idx), 'k-', 'LineWidth', 1.5); hold on;
    plot(time, Policy.tau_x_path(valid_idx), 'c--', 'LineWidth', 1.5);
    legend('\tau_{int}', '\tau_{x}');
    title('9. Other Wedges'); grid on; xlim([0 T]);
end