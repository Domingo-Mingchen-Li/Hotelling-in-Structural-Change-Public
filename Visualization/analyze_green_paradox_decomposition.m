function analyze_green_paradox_decomposition()
    % ANALYZE_GREEN_PARADOX_DECOMPOSITION_3WAY
    % 
    % e_t / e_{t+k} = [h_{t+k} / h_t] * [Y_t / Y_{t+k}] * [\bar{\beta}_t / \bar{\beta}_{t+k}]
    % \Delta \ln(e_1/e_{10}) = \Delta \ln(Price) + \Delta \ln(Macro_Growth) + \Delta \ln(Pure_Structural)

    clc; close all;
    fprintf('Initializing 3-Way Decomposition Engine...\n');

    % Input the path to any MATLAB data file (.mat) located in the Sim_Results directory.
    file_path = 'Sim_Results\Counterfactual_of_FrozenStructuralChange.mat';
    
    if ~exist(file_path, 'file')
        error('File not found. Please verify the data path.');
    end
    
    fprintf('Loading data...\n');
    load(file_path, 'SimData');
    
    Base = SimData.BaseRes;
    Final = SimData.FinalRes;
    p = SimData.p;
    Pol = SimData.TargetPol;
    BasePol = SimData.BasePol; 

    beta_m = p.beta_m;
    beta_s = p.beta_s;
    beta_x = p.beta_x;

    idx_1 = 1;
    idx_10 = 10;

    function [e_ratio, price_ratio, growth_ratio, struct_ratio] = get_components(Res, Policy)
        Rev_m = Res.p_m .* Res.c_m;
        Rev_s = Res.p_s .* Res.c_s;
        Rev_x = Res.I; 
        
        Y = Rev_m + Rev_s + Rev_x;
        
        E_tilde = beta_m * Rev_m + beta_s * Rev_s + beta_x * Rev_x;
        
        beta_bar = E_tilde ./ Y;

        growth_ratio = Y(idx_1) / Y(idx_10);
        struct_ratio = beta_bar(idx_1) / beta_bar(idx_10);

        tau_e_1 = Policy.tau_e_path(idx_1);
        tau_e_10 = Policy.tau_e_path(idx_10);
        h_gross_1 = Res.h(idx_1) * (1 + tau_e_1);
        h_gross_10 = Res.h(idx_10) * (1 + tau_e_10);
        price_ratio = h_gross_10 / h_gross_1;

        e_ratio = Res.R_input(idx_1) / Res.R_input(idx_10);
    end

    [e_rat_base, p_rat_base, g_rat_base, s_rat_base] = get_components(Base, BasePol);
    [e_rat_pol,  p_rat_pol,  g_rat_pol,  s_rat_pol]  = get_components(Final, Pol);

    log_e_diff = log(e_rat_pol) - log(e_rat_base);
    log_p_diff = log(p_rat_pol) - log(p_rat_base);
    log_g_diff = log(g_rat_pol) - log(g_rat_base);
    log_s_diff = log(s_rat_pol) - log(s_rat_base);
    
    residual = log_e_diff - (log_p_diff + log_g_diff + log_s_diff);

    fprintf('\n[Analysis Horizon]: t=0 to t=9\n');
    fprintf('------------------------------------------------------\n');
    fprintf('Total Tilt (Delta ln(e1/e10))  : %8.4f\n', log_e_diff);
    fprintf('   |= 1. Hotelling Effect      : %8.4f  (%5.1f%%)\n', log_p_diff, 100 * log_p_diff / log_e_diff);
    fprintf('   |= 2. Growth Effect         : %8.4f  (%5.1f%%)\n', log_g_diff, 100 * log_g_diff / log_e_diff);
    fprintf('   |= 3. Structural Effect     : %8.4f  (%5.1f%%)\n', log_s_diff, 100 * log_s_diff / log_e_diff);
    fprintf('   |= Math Residual            : %8.4e\n', residual);
    fprintf('------------------------------------------------------\n');

    figure('Name', '3-Way Decomposition of Green Paradox', 'Color', 'w', 'Position', [300, 200, 700, 500]);
    
    categories = categorical({'Total Tilt', 'Hotelling Effect', 'Growth Effect', 'Structural Change Effect'});
    categories = reordercats(categories, {'Total Tilt', 'Hotelling Effect', 'Growth Effect', 'Structural Change Effect'});
    values = [log_e_diff, log_p_diff, log_g_diff, log_s_diff] * 100; 
    
    b = bar(categories, values, 'FaceColor', 'flat');
    b.CData(1,:) = [0.2, 0.2, 0.5]; 
    b.CData(2,:) = [0.8, 0.3, 0.3]; 
    b.CData(3,:) = [0.9, 0.6, 0.2]; 
    b.CData(4,:) = [0.2, 0.6, 0.3]; 
    
    title('3-Way Decomposition of the Green Paradox (Log-Difference %)', 'FontSize', 14);
    ylabel('Relative Change vs Baseline (%)', 'FontSize', 12);
    grid on;
    
    for i = 1:numel(values)
        text(i, values(i) + sign(values(i))*0.5, sprintf('%.1f%%', values(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
    end
end