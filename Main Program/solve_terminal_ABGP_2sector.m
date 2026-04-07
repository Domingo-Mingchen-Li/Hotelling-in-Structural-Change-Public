function term = solve_terminal_ABGP_2sector(K_T, R_T, params, T_idx, Policy)
    get_pol = @(path) path(min(length(path), T_idx));
    
    TermPol.tau_c_m = get_pol(Policy.tau_c_m_path);
    TermPol.tau_c_s = get_pol(Policy.tau_c_s_path);
    TermPol.tau_x   = get_pol(Policy.tau_x_path);
    TermPol.tau_int = get_pol(Policy.tau_int_path);
    TermPol.tau_e   = get_pol(Policy.tau_e_path);
    
    g = get_growth_rates_2sector(params); 
    
    A_T.m = params.A0.m * (params.g_A_m ^ T_idx);
    A_T.s = params.A0.s * (params.g_A_s ^ T_idx);
    A_T.x = params.A0.x * (params.g_A_x ^ T_idx);
    A_T.e = params.A0.e * (params.g_A_e ^ T_idx);
    
    if g.g_e >= 1
        fprintf('Warning: g_e >= 1 (%.4f), Resource sum diverges.\n', g.g_e);
    end
    e_T = R_T * (1 - g.g_e);
    
    net_growth_K = g.g_K - (1 - params.delta);
    I_T = net_growth_K * K_T;
    
    if I_T < 0
        fprintf('Warning: Terminal Investment is negative (g_K=%.4f).\n', g.g_K);
    end
    
    r_target = g.r_ss;

    obj_fun = @(E_guess) check_r_target(E_guess, K_T, e_T, I_T, A_T, r_target, params, TermPol);
    
    E_min = 1e-3 * K_T;
    E_max = 100.0 * K_T; 
    
    options = optimset('Display','off', 'TolX', 1e-14, 'TolFun', 1e-14);
    
    try
        [E_T, ~, exitflag] = fzero(obj_fun, [E_min, E_max], options);
    catch
        options_min = optimset('Display','off', 'TolX', 1e-14, 'TolFun', 1e-14);
        obj_sq = @(E) (check_r_target(E, K_T, e_T, I_T, A_T, r_target, params, TermPol))^2;
        E_T = fminbnd(obj_sq, E_min, E_max, options_min);
        exitflag = 1; 
    end

    if exitflag < 0
        fprintf('Warning: ABGP Solver convergence failed. Falling back to bounded result.\n');
    end

    [stat] = static_kernel_2sector(K_T, e_T, E_T, I_T, A_T, params, TermPol);
    
    term.K = K_T;
    term.R = R_T;
    term.E = E_T;
    term.e = e_T; 
    term.I = I_T; 
    
    term.r = stat.r;
    term.w = stat.w;
    term.h = stat.h; 
    
    term.p_int = stat.p_int;
    term.p_m   = stat.p_m;
    term.p_s   = stat.p_s;
    
    term.P_star = stat.P_star; 
    
    term.Y_GDP = E_T + I_T * (1 + TermPol.tau_x);
end

function diff = check_r_target(E_guess, K_val, e_val, I_phys, A_val, r_target, p, pol)
    if E_guess <= 0, diff = 1e10; return; end
    [stat] = static_kernel_2sector(K_val, e_val, E_guess, I_phys, A_val, p, pol);
    diff = stat.r - r_target;
end