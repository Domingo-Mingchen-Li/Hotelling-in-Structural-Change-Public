function g = get_growth_rates_2sector(p)
    be = p.sect.e.beta; ae = p.sect.e.alpha;
    bx = p.sect.x.beta; ax = p.sect.x.alpha;
    eps_val = p.epsilon;
    
    w_alpha = p.omega_m*p.sect.m.alpha + p.omega_s*p.sect.s.alpha;
    w_beta  = p.omega_m*p.sect.m.beta  + p.omega_s*p.sect.s.beta;

    term_supply = (1/be) - (ae * bx) / (be * (1-ax));
    term_demand = eps_val * ( (bx / (1-ax)) * w_alpha + w_beta );
    phi = term_supply - term_demand;
    
    exp_beta = 1 / phi;
    exp_Ae   = 1 / (be * phi);
    exp_Ax   = (eps_val / (1-ax)) / phi; 
    
    exp_Am   = (p.omega_m * eps_val) / phi;
    exp_As   = (p.omega_s * eps_val) / phi;
    
    g_pivot = (p.beta ^ exp_beta) * ...
              (p.g_A_e ^ exp_Ae) * ...
              (p.g_A_x ^ exp_Ax) * ...
              (p.g_A_m ^ exp_Am) * (p.g_A_s ^ exp_As);

    g_K = (p.g_A_x ^ (1/(1-ax))) * (g_pivot ^ (bx/(1-ax)));
    g_E = g_K;
    
    exp_e_Ax = - ae / ((1-ax) * be); 
    exp_e_pivot = (1/be) * (1 - (bx * ae)/(1-ax));
    
    g_e = (p.g_A_e ^ (-1/be)) * ...
          (p.g_A_x ^ exp_e_Ax) * ...
          (g_pivot ^ exp_e_pivot);
          
    exp_h_Ax = (1/(1-ax)) * (1 + ae/be);
    exp_h_pivot = (bx/(1-ax)) * (1 + ae/be) - (1/be);
    
    g_h = (p.g_A_e ^ (1/be)) * ...
          (p.g_A_x ^ exp_h_Ax) * ...
          (g_pivot ^ exp_h_pivot);

    calc_gpi = @(gAi, ai, bi) (gAi^-1) * (p.g_A_x^((1-ai)/(1-ax))) * (g_pivot ^ (bx*(1-ai)/(1-ax) - bi));
    
    g.g_pm = calc_gpi(p.g_A_m, p.sect.m.alpha, p.sect.m.beta);
    g.g_ps = calc_gpi(p.g_A_s, p.sect.s.alpha, p.sect.s.beta);
    
    g.g_P = (g.g_pm^p.omega_m) * (g.g_ps^p.omega_s);

    R_gross = (1/p.beta) * (g_E^(1-eps_val)) * (g.g_P^eps_val);
    g.r_ss = R_gross - (1 - p.delta);
    
    g.g_m = g_pivot; 
    g.g_K = g_K; 
    g.g_E = g_E; 
    g.g_e = g_e; 
    g.g_h = g_h;
end