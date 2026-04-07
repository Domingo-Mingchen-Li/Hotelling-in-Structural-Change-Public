function test_counterfactual()
    clc; clear; close all;
    fprintf('Initializing Counterfactual Test: Frozen Structural Change\n');

    p = struct();
    p.epsilon = 0.22; p.beta = 0.96; p.delta = 0.05;   
    
    p.alpha_m = 0.4146; p.beta_m = 0.3024; 
    p.alpha_s = 0.4452; p.beta_s = 0.1035; 
    p.alpha_x = 0.2223; p.beta_x = 0.5415; 
    p.alpha_e = 0.3506; p.beta_e = 0.4737;
    
    p.g_A_m = 1.0331; p.g_A_s = 1.0000; 
    p.g_A_x = 1.0352; p.g_A_e = 1.0118;

    make_sect = @(a, b) struct('alpha', a, 'beta', b);
    p.sect.m = make_sect(p.alpha_m, p.beta_m);
    p.sect.s = make_sect(p.alpha_s, p.beta_s);
    p.sect.x = make_sect(p.alpha_x, p.beta_x);
    p.sect.e = make_sect(p.alpha_e, p.beta_e);
    
    p.A0.m = 1.0; p.A0.s = 1.0; p.A0.x = 1.0; p.A0.e = 1.0;

    p.omega_m = 0.10; 
    p.omega_s = 0.90; 
    
    p.xi = 1;      
    p.zeta_m = 0.29807;   
    p.zeta_s = -0.29807;  
    
    T_periods = 75;

    Targets.K0 = 5.8439; 
    Targets.R0 = 1137.4505; 
    
    ShootingSettings = struct();
    ShootingSettings.T_periods  = T_periods;
    ShootingSettings.tol        = 1e-4;    
    ShootingSettings.max_iter   = 50;     
    ShootingSettings.perturb_K  = 1e-9;   
    ShootingSettings.perturb_R  = 1e-12;   
    ShootingSettings.min_lambda = 1e-2;   

    fprintf('Establishing baseline...\n');
    
    BasePol = struct();
    zeros_path = zeros(1, T_periods + 1);
    BasePol.tau_c_m_path = zeros_path;
    BasePol.tau_c_s_path = zeros_path;
    BasePol.tau_x_path   = zeros_path;
    BasePol.tau_int_path = zeros_path;
    BasePol.tau_e_path   = zeros_path;
    
    Guess.KT = 1245.90562; 
    Guess.RT = 52.01390;
    
    try
        [BaseRes, base_KT, base_RT, success] = solve_shooting_2sector_inverted(p, BasePol, Guess, Targets, ShootingSettings);
    catch ME
        error('Baseline solver crashed: %s\n', ME.message);
    end
    
    if ~success
        error('Baseline establishment failed.');
    end
    fprintf('Baseline ready: KT=%.5f, RT=%.5f\n', base_KT, base_RT);
    
    BaselineStart.KT = base_KT;
    BaselineStart.RT = base_RT;
    BaselineStart.K  = BaseRes.K;
    BaselineStart.R  = BaseRes.R;

    val_m_0 = BaseRes.p_m(1) * (1 + BasePol.tau_c_m_path(1)) * BaseRes.c_m(1);
    val_s_0 = BaseRes.p_s(1) * (1 + BasePol.tau_c_s_path(1)) * BaseRes.c_s(1);
    total_val_0 = val_m_0 + val_s_0;
    
    real_share_m_0 = val_m_0 / total_val_0;
    real_share_s_0 = val_s_0 / total_val_0;
    
    fprintf('Extracted baseline expenditure shares: omega_m = %.4f, omega_s = %.4f\n', real_share_m_0, real_share_s_0);

    fprintf('Configuring counterfactual targets...\n');
    
    TargetParams = struct();
    
    TargetParams.tau_e_0 = 0.00;  
    TargetParams.g_tau_e = 0.00;  
    TargetParams.tau_int = 0.00;  
    TargetParams.tau_x   = 0.00;  
    TargetParams.tau_c_m = 0.00;  
    TargetParams.tau_c_s = 0.00;  
    
    % To freeze structural change completely: 
    % 1. Set TargetParams.zeta_m = 0
    % 2. Uncomment the following two lines for omega_m and omega_s
    TargetParams.zeta_m  = 0.8*p.zeta_m;  
    %TargetParams.omega_m = real_share_m_0;
    %TargetParams.omega_s = real_share_s_0;
    
    TargetParams.omega_m = p.omega_m;
    TargetParams.omega_s = p.omega_s;
    
    fprintf('Target Policies: zeta_m -> %.4f, omega_m -> %.4f\n', TargetParams.zeta_m, TargetParams.omega_m);

    fprintf('Initiating precision homotopy...\n');
    
    [FinalRes, History] = run_precision_homotopy_2sector(p, Targets, BaselineStart, TargetParams, ShootingSettings);
    
    if ~isempty(FinalRes)
        TargetPol = reconstruct_target_policy(TargetParams, T_periods);
        
        visualize_all(BaseRes, FinalRes, BasePol, TargetPol, T_periods);
        
        save_dir = 'Sim_Results';
        if ~exist(save_dir, 'dir')
            mkdir(save_dir);
        end
        
        file_name = 'Counterfactual_FrozenStructChange.mat';
                        
        full_path = fullfile(save_dir, file_name);
        
        SimData = struct();
        SimData.p = p;
        SimData.Targets = Targets;
        SimData.BasePol = BasePol;
        SimData.TargetPol = TargetPol;
        SimData.TargetParams = TargetParams;
        SimData.BaseRes = BaseRes;
        SimData.FinalRes = FinalRes;
        SimData.History = History;
        SimData.T_periods = T_periods;
        
        save(full_path, 'SimData');
        fprintf('Results saved to: %s\n', full_path);
        
    else
        fprintf('Homotopy incomplete. Plotting and saving aborted.\n');
    end
end

function Pol = reconstruct_target_policy(Target, T)
    time_vec = 0:T;
    tau_0 = 0; g = 0;
    if isfield(Target, 'tau_e_0'), tau_0 = Target.tau_e_0; end
    if isfield(Target, 'g_tau_e'), g = Target.g_tau_e; end
    Pol.tau_e_path = tau_0 * ((1 + g) .^ time_vec);
    
    simple_taxes = {'tau_x', 'tau_int', 'tau_c_m', 'tau_c_s'};
    for i = 1:length(simple_taxes)
        name = simple_taxes{i};
        path_name = [name '_path'];
        val = 0;
        if isfield(Target, name), val = Target.(name); end
        Pol.(path_name) = val * ones(1, T+1);
    end
end

function visualize_all(Base, Final, BasePol, TargetPol, T)
    obs_period = 40; 
    time = 0:T;
    

    Val_m_base = Base.p_m .* (1 + BasePol.tau_c_m_path) .* Base.c_m;
    Val_s_base = Base.p_s .* (1 + BasePol.tau_c_s_path) .* Base.c_s;
    C_in_base = Val_m_base + Val_s_base;
    Share_Cm_base = Val_m_base ./ C_in_base;
    Share_Cs_base = Val_s_base ./ C_in_base;
    
    Val_m_pol = Final.p_m .* (1 + TargetPol.tau_c_m_path) .* Final.c_m;
    Val_s_pol = Final.p_s .* (1 + TargetPol.tau_c_s_path) .* Final.c_s;
    C_in_pol = Val_m_pol + Val_s_pol;
    Share_Cm_pol = Val_m_pol ./ C_in_pol;
    Share_Cs_pol = Val_s_pol ./ C_in_pol;
    
    L_in_base = Base.l_m + Base.l_s;
    Share_Lm_base = Base.l_m ./ L_in_base;
    Share_Ls_base = Base.l_s ./ L_in_base;
    
    L_in_pol = Final.l_m + Final.l_s;
    Share_Lm_pol = Final.l_m ./ L_in_pol;
    Share_Ls_pol = Final.l_s ./ L_in_pol;

    Inv_exp_base = Base.I .* (1 + BasePol.tau_x_path);
    Inv_exp_final = Final.I .* (1 + TargetPol.tau_x_path);
    Inv_rate_base = Inv_exp_base ./ Base.Y_GDP;
    Inv_rate_pol = Inv_exp_final ./ Final.Y_GDP;
    
    Exp_share_base = Base.E ./ Base.Y_GDP;
    Exp_share_pol = Final.E ./ Final.Y_GDP;


    figure('Name', 'Figure 1: Structural Transformation (Absolute)', 'Color', 'w', 'Position', [50, 500, 1200, 600]);
    subplot(2,2,1); plot(time, Share_Cm_base, 'k--', 'LineWidth', 1.5); hold on; plot(time, Share_Cm_pol, 'b-', 'LineWidth', 2); title('Consumption Share: Manufacturing'); legend('Baseline', 'Frozen \zeta & \omega'); grid on; ylim([0 1]);
    subplot(2,2,2); plot(time, Share_Cs_base, 'k--', 'LineWidth', 1.5); hold on; plot(time, Share_Cs_pol, 'r-', 'LineWidth', 2); title('Consumption Share: Services'); legend('Baseline', 'Frozen \zeta & \omega'); grid on; ylim([0 1]);
    subplot(2,2,3); plot(time, Share_Lm_base, 'k--', 'LineWidth', 1.5); hold on; plot(time, Share_Lm_pol, 'b-', 'LineWidth', 2); title('Labor Share: Manuf (Lm/L_{cons})'); grid on; ylim([0 1]);
    subplot(2,2,4); plot(time, Share_Ls_base, 'k--', 'LineWidth', 1.5); hold on; plot(time, Share_Ls_pol, 'r-', 'LineWidth', 2); title('Labor Share: Services (Ls/L_{cons})'); grid on; ylim([0 1]);
    set(findobj(gcf, 'Type', 'axes'), 'XLim', [0, obs_period]);

    figure('Name', 'Figure 2: Macro Dynamics & Tech (Absolute)', 'Color', 'w', 'Position', [100, 400, 1400, 700]);
    subplot(2,3,1); plot(time, Base.I, 'k--', 'LineWidth', 1.5); hold on; plot(time, Final.I, 'g-', 'LineWidth', 2); title('Investment Level (I)'); legend('Baseline', 'Frozen Structure'); grid on;
    subplot(2,3,2); plot(time, Base.E, 'k--', 'LineWidth', 1.5); hold on; plot(time, Final.E, 'b-', 'LineWidth', 2); title('Nominal Expenditure (E)'); legend('Baseline', 'Frozen Structure'); grid on;
    subplot(2,3,3); plot(time, Base.p_m ./ Base.p_s, 'k--', 'LineWidth', 1.5); hold on; plot(time, Final.p_m ./ Final.p_s, 'm-', 'LineWidth', 2); title('Relative Price (Pm / Ps)'); legend('Baseline', 'Frozen Structure'); grid on;
    subplot(2,3,4); plot(time, Inv_rate_base, 'k--', 'LineWidth', 1.5); hold on; plot(time, Inv_rate_pol, 'g-', 'LineWidth', 2); title('Inv. Rate (I / GDP)'); grid on;
    subplot(2,3,5); plot(time, Exp_share_base, 'k--', 'LineWidth', 1.5); hold on; plot(time, Exp_share_pol, 'b-', 'LineWidth', 2); title('Exp. Share (E / GDP)'); grid on;
    subplot(2,3,6); plot(time, Final.A_m, 'b-o', 'MarkerIndices', 1:10:T, 'MarkerSize', 4); hold on; plot(time, Final.A_s, 'r-'); plot(time, Final.A_x, 'g-'); plot(time, Final.A_e, 'm-'); title('Tech Levels (A)'); legend('Am','As','Ax','Ae'); grid on; set(gca,'YScale','log');
    set(findobj(gcf, 'Type', 'axes'), 'XLim', [0, obs_period]);

    figure('Name', 'Figure 3: State Variables & Extraction (Absolute)', 'Color', 'w', 'Position', [150, 100, 1200, 400]);
    subplot(1,3,1); plot(time, Base.K, 'k--', 'LineWidth', 1.5); hold on; plot(time, Final.K, 'b-', 'LineWidth', 2); title('Capital Stock (K)'); legend('Baseline', 'Frozen Structure'); grid on; xlabel('Time');
    subplot(1,3,2); plot(time, Base.R, 'k--', 'LineWidth', 1.5); hold on; plot(time, Final.R, 'r-', 'LineWidth', 2); title('Resource Stock (R)'); legend('Baseline', 'Frozen Structure'); grid on; xlabel('Time');
    subplot(1,3,3); plot(time, Base.R_input, 'k--', 'LineWidth', 1.5); hold on; plot(time, Final.R_input, 'm-', 'LineWidth', 2); title('Extraction Flow (R_{input})'); legend('Baseline', 'Frozen Structure'); grid on; xlabel('Time');
    set(findobj(gcf, 'Type', 'axes'), 'XLim', [0, obs_period]);

    figure('Name', 'Figure 4: Structural Transformation (Deviation from Baseline)', 'Color', 'w', 'Position', [200, 500, 1200, 600]);
    subplot(2,2,1); plot(time, (Share_Cm_pol - Share_Cm_base) * 100, 'b-', 'LineWidth', 2); yline(0, 'k--'); title('\Delta Cons Share: Manuf (pp)'); grid on; ylabel('Percentage Points');
    subplot(2,2,2); plot(time, (Share_Cs_pol - Share_Cs_base) * 100, 'r-', 'LineWidth', 2); yline(0, 'k--'); title('\Delta Cons Share: Services (pp)'); grid on; ylabel('Percentage Points');
    subplot(2,2,3); plot(time, (Share_Lm_pol - Share_Lm_base) * 100, 'b-', 'LineWidth', 2); yline(0, 'k--'); title('\Delta Labor Share: Manuf (pp)'); grid on; ylabel('Percentage Points');
    subplot(2,2,4); plot(time, (Share_Ls_pol - Share_Ls_base) * 100, 'r-', 'LineWidth', 2); yline(0, 'k--'); title('\Delta Labor Share: Services (pp)'); grid on; ylabel('Percentage Points');
    set(findobj(gcf, 'Type', 'axes'), 'XLim', [0, obs_period]);

    figure('Name', 'Figure 5: Macro Dynamics (Deviation from Baseline)', 'Color', 'w', 'Position', [250, 400, 1400, 400]);
    subplot(1,4,1); plot(time, (Final.I ./ Base.I - 1) * 100, 'g-', 'LineWidth', 2); yline(0, 'k--'); title('\Delta Investment Level (%)'); grid on; ylabel('% Deviation');
    subplot(1,4,2); plot(time, (Final.Y_GDP ./ Base.Y_GDP - 1) * 100, 'b-', 'LineWidth', 2); yline(0, 'k--'); title('\Delta GDP Level (%)'); grid on; ylabel('% Deviation');
    subplot(1,4,3); plot(time, ((Final.p_m ./ Final.p_s) ./ (Base.p_m ./ Base.p_s) - 1) * 100, 'm-', 'LineWidth', 2); yline(0, 'k--'); title('\Delta Relative Price Pm/Ps (%)'); grid on; ylabel('% Deviation');
    subplot(1,4,4); plot(time, (Inv_rate_pol - Inv_rate_base) * 100, 'g-', 'LineWidth', 2); yline(0, 'k--'); title('\Delta Inv. Rate (I/GDP) (pp)'); grid on; ylabel('Percentage Points');
    set(findobj(gcf, 'Type', 'axes'), 'XLim', [0, obs_period]);

    figure('Name', 'Figure 6: State Variables & Extraction (Deviation from Baseline)', 'Color', 'w', 'Position', [300, 100, 1200, 400]);
    subplot(1,3,1); plot(time, (Final.K ./ Base.K - 1) * 100, 'b-', 'LineWidth', 2); yline(0, 'k--'); title('\Delta Capital Stock K (%)'); grid on; ylabel('% Deviation');
    subplot(1,3,2); plot(time, (Final.R ./ Base.R - 1) * 100, 'r-', 'LineWidth', 2); yline(0, 'k--'); title('\Delta Resource Stock R (%)'); grid on; ylabel('% Deviation');
    subplot(1,3,3); plot(time, (Final.R_input ./ Base.R_input - 1) * 100, 'm-', 'LineWidth', 2); yline(0, 'k--'); title('\Delta Extraction Flow R_{input} (%)'); grid on; ylabel('% Deviation');
    set(findobj(gcf, 'Type', 'axes'), 'XLim', [0, obs_period]);
end