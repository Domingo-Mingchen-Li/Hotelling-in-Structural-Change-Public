function [solution, success] = solve_backward_step_2sector(future_state, params, tech, current_guesses, Policy, t_idx)
    if nargin < 5, Policy = struct(); end
    if ~isfield(Policy, 'tau_c_m'), Policy.tau_c_m = 0; end
    if ~isfield(Policy, 'tau_c_s'), Policy.tau_c_s = 0; end
    if ~isfield(Policy, 'tau_x'),   Policy.tau_x = 0; end
    if ~isfield(Policy, 'tau_int'), Policy.tau_int = 0; end
    if ~isfield(Policy, 'tau_e'),   Policy.tau_e = 0; end
    if ~isfield(Policy, 'tau_x_next'), Policy.tau_x_next = 0; end

    h_net_current = future_state.h_next / (1 - params.delta + future_state.r_next);
    
    if nargin < 4 || isempty(current_guesses), r0 = 0.05; else, r0 = current_guesses.r; end
    
    r_min = 1e-5; r_max = 5.0; 
    options = optimset('TolX', 1e-8, 'Display', 'off');
    
    obj_fun = @(r) check_capital_market(r, h_net_current, future_state, params, tech, Policy);
    
    try
        [r_sol, ~, exitflag] = fzero(obj_fun, r0, options);
        if exitflag > 0
            success = true;
        else
            if obj_fun(r_min) * obj_fun(r_max) < 0
                 r_sol = fzero(obj_fun, [r_min, r_max], options); success = true;
            else
                 r_sol = r0; success = false;
            end
        end
    catch
        r_sol = r0; success = false;
    end
    
    if success
        solution = compute_full_equilibrium(r_sol, h_net_current, future_state, params, tech, Policy);
    else
        solution = struct(); 
        try
            solution = compute_full_equilibrium(r_sol, h_net_current, future_state, params, tech, Policy);
        catch
        end
    end
end

function diff_K = check_capital_market(r_guess, h_net, future, params, tech, Pol)
    if r_guess <= 0, diff_K = 1e10; return; end

    [w, f_prod, f_user, p_m, p_s] = reconstruct_prices(r_guess, h_net, params, tech, Pol);
    
    if ~isreal(w) || w <= 0 || f_prod <= 0, diff_K = 1e10; return; end

    p_m_gross = p_m * (1 + Pol.tau_c_m);
    p_s_gross = p_s * (1 + Pol.tau_c_s);
    
    P_star = (p_m_gross ^ params.omega_m) * (p_s_gross ^ params.omega_s);
    
    term_future = (future.E_next^(params.epsilon - 1)) / ...
                  (future.P_star_next^params.epsilon * (1 + Pol.tau_x_next));
              
    euler_rhs = params.beta * (1 - params.delta + future.r_next) * term_future;
    
    lhs_val = euler_rhs * (P_star^params.epsilon) * (1 + Pol.tau_x);
    E_current = lhs_val ^ (1 / (params.epsilon - 1));

    [c_m, c_s] = compute_consumption(E_current, p_m_gross, p_s_gross, P_star, params);

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
    
    l_x = (1 - L_cons - Const_LE) / (1 + Coeff_LE_LX);
    
    if l_x < 0, diff_K = 1e5 + abs(l_x)*1e5; return; end
    
    I_real = w * l_x / gamma_x; 
    
    K_supply = (future.K_next - I_real) / (1 - params.delta);
    if K_supply <= 0, diff_K = 1e8; return; end
    
    val_K_m = params.alpha_m * p_m * c_m;
    val_K_s = params.alpha_s * p_s * c_s;
    val_K_x = params.alpha_x * 1 * I_real;
    
    m_x_real = params.beta_x * I_real / f_user;
    Y_int_total = M_cons + m_x_real;
    val_K_e = params.alpha_e * f_prod * Y_int_total; 
    
    K_demand = (val_K_m + val_K_s + val_K_x + val_K_e) / r_guess;
    
    diff_K = K_supply - K_demand;
end

function [w, f_prod, f_user, p_m, p_s] = reconstruct_prices(r, h_net, p, tech, Pol)
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
    
    ln_w = -(term_const + term_r) / coeff_w;
    w = exp(ln_w);
    
    ln_f_prod = C_e + p.alpha_e*log(r) + gamma_e*ln_w + p.beta_e*cost_h_log;
    f_prod = exp(ln_f_prod);
    f_user = f_prod * (1 + Pol.tau_int);
    
    calc_p = @(C, al, be, ga) exp(C + al*log(r) + ga*log(w) + be*log(f_user));
    
    C_m = get_log_C(tech.A_m, p.alpha_m, p.beta_m, 1-p.alpha_m-p.beta_m);
    p_m = calc_p(C_m, p.alpha_m, p.beta_m, 1-p.alpha_m-p.beta_m);
    
    C_s = get_log_C(tech.A_s, p.alpha_s, p.beta_s, 1-p.alpha_s-p.beta_s);
    p_s = calc_p(C_s, p.alpha_s, p.beta_s, 1-p.alpha_s-p.beta_s);
end

function [c_m, c_s] = compute_consumption(E, p_m_g, p_s_g, P_star, p)
    exp_m = p.epsilon * p.omega_m + p.zeta_m; 
    exp_s = p.epsilon * p.omega_s + p.zeta_s;
    
    Pi_star = (p_m_g ^ exp_m) * (p_s_g ^ exp_s);
    
    term_pigl = p.xi * (E^(1 - p.epsilon)) * Pi_star * (P_star^p.epsilon);
    
    val_m = p.omega_m * E + p.zeta_m * term_pigl;
    val_s = p.omega_s * E + p.zeta_s * term_pigl;
    
    c_m = val_m / p_m_g;
    c_s = val_s / p_s_g;
    
    c_m = max(c_m, 1e-9); c_s = max(c_s, 1e-9);
end

function sol = compute_full_equilibrium(r, h_net, future, params, tech, Pol)
    [w, f_prod, f_user, p_m, p_s] = reconstruct_prices(r, h_net, params, tech, Pol);
    
    p_m_gross = p_m * (1 + Pol.tau_c_m);
    p_s_gross = p_s * (1 + Pol.tau_c_s);
    P_star = (p_m_gross ^ params.omega_m) * (p_s_gross ^ params.omega_s);
    
    term_future = (future.E_next^(params.epsilon - 1)) / (future.P_star_next^params.epsilon * (1 + Pol.tau_x_next));
    euler_rhs = params.beta * (1 - params.delta + future.r_next) * term_future;
    lhs_val = euler_rhs * (P_star^params.epsilon) * (1 + Pol.tau_x);
    E = lhs_val ^ (1 / (params.epsilon - 1));
    
    [c_m, c_s] = compute_consumption(E, p_m_gross, p_s_gross, P_star, params);
    
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
    l_x = (1 - L_cons - Const_LE) / (1 + Coeff_LE_LX);
    
    I_real = w * l_x / gamma_x;
    
    m_x_real = params.beta_x * I_real / f_user;
    Y_int_total = M_cons + m_x_real;
    
    l_e = gamma_e * f_prod * Y_int_total / w;
    
    h_gross = h_net * (1 + Pol.tau_e);
    R_input = params.beta_e * f_prod * Y_int_total / h_gross;
    
    K_curr = (future.K_next - I_real) / (1 - params.delta);
    R_stock_curr = future.R_next + R_input;
    
    sol.r = r; sol.w = w; sol.h = h_net; 
    sol.p_int = f_prod; sol.p_int_user = f_user; 
    sol.p_m = p_m; sol.p_s = p_s;
    sol.E = E; sol.I = I_real; 
    sol.R_input = R_input; 
    sol.K = K_curr; sol.R = R_stock_curr;
    sol.c_m = c_m; sol.c_s = c_s;
    sol.l_m = val_L_m/w; sol.l_s = val_L_s/w; sol.l_x = l_x; sol.l_e = l_e;
    sol.P_star = P_star; 
end