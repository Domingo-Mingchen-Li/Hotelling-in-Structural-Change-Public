function [stat] = static_kernel_2sector(K_supply, e_supply, E_nom, I_phys, A_tech, params, Pol)
    Rev_m = (params.omega_m * E_nom) / (1 + Pol.tau_c_m);
    Rev_s = (params.omega_s * E_nom) / (1 + Pol.tau_c_s);
    Rev_x = I_phys; 
    
    Total_Rev = Rev_m + Rev_s + Rev_x;
    
    w_m = Rev_m / Total_Rev;
    w_s = Rev_s / Total_Rev;
    w_x = Rev_x / Total_Rev;
    
    pass_through = 1 / (1 + Pol.tau_int);
    
    sm = params.sect.m; ss = params.sect.s; 
    sx = params.sect.x; se = params.sect.e; 
    
    get_shares = @(s) struct(...
        'K', s.alpha + s.beta * pass_through * se.alpha, ...
        'L', (1 - s.alpha - s.beta) + s.beta * pass_through * (1 - se.alpha - se.beta), ...
        'R', s.beta * pass_through * se.beta ...
    );

    PSI.m = get_shares(sm);
    PSI.s = get_shares(ss);
    PSI.x = get_shares(sx);
    
    Agg_K = w_m*PSI.m.K + w_s*PSI.s.K + w_x*PSI.x.K;
    Agg_L = w_m*PSI.m.L + w_s*PSI.s.L + w_x*PSI.x.L;
    Agg_R = w_m*PSI.m.R + w_s*PSI.s.R + w_x*PSI.x.R;
    
    Lambda_w = (Agg_L / Agg_K) * (K_supply / 1.0); 
    
    Lambda_h_gross = (Agg_R / Agg_K) * (K_supply / e_supply);
    
    get_log_C = @(A, al, be, ga) -log(A) - al*log(al) - be*log(be) - ga*log(ga);
    
    C_x = get_log_C(A_tech.x, sx.alpha, sx.beta, 1-sx.alpha-sx.beta);
    C_e = get_log_C(A_tech.e, se.alpha, se.beta, 1-se.alpha-se.beta);
    
    lx_x = 1 - sx.alpha - sx.beta;
    le_e = 1 - se.alpha - se.beta;
    
    ln_lambda_p_int = C_e + le_e*log(Lambda_w) + se.beta*log(Lambda_h_gross);
    
    wedge_int_log = log(1 + Pol.tau_int);
    ln_lambda_p_int_user = ln_lambda_p_int + wedge_int_log;
    
    numer = C_x + lx_x*log(Lambda_w) + sx.beta*ln_lambda_p_int_user;
    ln_r = -numer;
    r = exp(ln_r);
    
    w = r * Lambda_w;
    h_gross = r * Lambda_h_gross;
    h_net = h_gross / (1 + Pol.tau_e); 
    
    p_int = r * exp(ln_lambda_p_int);
    p_int_user = p_int * (1 + Pol.tau_int);
    
    calc_p = @(C, al, be, ga) exp(C + al*log(r) + ga*log(w) + be*log(p_int_user));
    
    C_m = get_log_C(A_tech.m, sm.alpha, sm.beta, 1-sm.alpha-sm.beta);
    p_m = calc_p(C_m, sm.alpha, sm.beta, 1-sm.alpha-sm.beta);
    
    C_s = get_log_C(A_tech.s, ss.alpha, ss.beta, 1-ss.alpha-ss.beta);
    p_s = calc_p(C_s, ss.alpha, ss.beta, 1-ss.alpha-ss.beta);
    
    p_m_g = p_m * (1 + Pol.tau_c_m);
    p_s_g = p_s * (1 + Pol.tau_c_s);
    P_star = (p_m_g ^ params.omega_m) * (p_s_g ^ params.omega_s);
    
    stat.r = r; stat.w = w; stat.h = h_net;
    stat.p_int = p_int; 
    stat.p_m = p_m; stat.p_s = p_s;
    stat.P_star = P_star;
    stat.Total_Rev = Total_Rev;
end