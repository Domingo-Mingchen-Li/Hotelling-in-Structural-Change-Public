function [FinalRes, History] = run_precision_homotopy_2sector(p, Targets, BaselineStart, TargetParams, ShootingSettings)
    clc;
    fprintf('Initializing Adaptive Homotopy Solver\n');

    if nargin < 5 || isempty(ShootingSettings)
        ShootingSettings = struct('T_periods', 100, 'tol', 1e-4, 'max_iter', 30, ...
                                  'perturb_K', 1e-4, 'perturb_R', 1e-5, 'min_lambda', 1e-8);
    end
    
    if isstruct(BaselineStart) && isfield(BaselineStart, 'KT')
        Current_KT = BaselineStart.KT; Current_RT = BaselineStart.RT;
    elseif isstruct(BaselineStart) && isfield(BaselineStart, 'K')
        Current_KT = BaselineStart.K(end); Current_RT = BaselineStart.R(end);
    else
        error('BaselineStart invalid.');
    end
    
    fprintf('>>> Solver: T=%d, MaxIter=%d, Perturb=%.1e/%.1e\n', ...
        ShootingSettings.T_periods, ShootingSettings.max_iter, ...
        ShootingSettings.perturb_K, ShootingSettings.perturb_R);
    
    base_zeta_m = p.zeta_m;
    if isfield(TargetParams, 'zeta_m')
        target_zeta_m = TargetParams.zeta_m;
        do_zeta_homotopy = true;
        fprintf('Homotopy enabled: zeta_m [%.4f -> %.4f]\n', base_zeta_m, target_zeta_m);
    else
        target_zeta_m = base_zeta_m;
        do_zeta_homotopy = false;
    end
    
    base_omega_m = p.omega_m;
    if isfield(TargetParams, 'omega_m')
        target_omega_m = TargetParams.omega_m;
        do_omega_homotopy = true;
        fprintf('Homotopy enabled: omega_m [%.4f -> %.4f]\n', base_omega_m, target_omega_m);
    else
        target_omega_m = base_omega_m;
        do_omega_homotopy = false;
    end
    
    base_R0 = Targets.R0;
    if isfield(TargetParams, 'target_R0')
        target_R0_val = TargetParams.target_R0;
        do_R0_homotopy = true;
        fprintf('Homotopy enabled: Target R0 [%.4f -> %.4f]\n', base_R0, target_R0_val);
    else
        target_R0_val = base_R0;
        do_R0_homotopy = false;
    end

    sectors = {'m', 's', 'x', 'e'};
    do_tech_homotopy = false;
    base_tech = struct();
    target_tech = struct();
    
    for i = 1:length(sectors)
        sec = sectors{i};
        
        base_tech.A0.(sec) = p.A0.(sec);
        field_A0 = ['target_A0_' sec];
        if isfield(TargetParams, field_A0)
            target_tech.A0.(sec) = TargetParams.(field_A0);
            do_tech_homotopy = true;
            fprintf('Homotopy enabled: A0_%s -> %.4f\n', sec, target_tech.A0.(sec));
        else
            target_tech.A0.(sec) = base_tech.A0.(sec);
        end
        
        base_g_name = ['g_A_' sec];
        base_tech.g_A.(sec) = p.(base_g_name);
        field_g_A = ['target_g_A_' sec];
        if isfield(TargetParams, field_g_A)
            target_tech.g_A.(sec) = TargetParams.(field_g_A);
            do_tech_homotopy = true;
            fprintf('Homotopy enabled: g_A_%s -> %.4f\n', sec, target_tech.g_A.(sec));
        else
            target_tech.g_A.(sec) = base_tech.g_A.(sec);
        end
    end
    
    current_progress = 0.0; 
    step_size = 0.005;       
    max_step_size = 0.005;   
    min_step_size = 1e-7;    
    
    History.Progress = 0;
    History.KT = Current_KT;
    History.RT = Current_RT;
    
    iter_count = 0; FinalRes = [];
    
    fprintf('Starting iteration: KT=%.4f, RT=%.4f\n', Current_KT, Current_RT);
    fprintf('Step size: %.2f%% (Max: %.2f%%)\n\n', step_size*100, max_step_size*100);
    
    while current_progress < 1.0
        iter_count = iter_count + 1;
        
        attempt_progress = current_progress + step_size;
        if attempt_progress > 1.0, attempt_progress = 1.0; end
        
        CurrentPolicy = generate_policy_from_params(TargetParams, attempt_progress, ShootingSettings.T_periods);
        
        p_curr = p;
        if do_zeta_homotopy
            curr_zeta_m = (1 - attempt_progress) * base_zeta_m + attempt_progress * target_zeta_m;
            p_curr.zeta_m = curr_zeta_m;
            p_curr.zeta_s = -curr_zeta_m; 
        end
        if do_omega_homotopy
            curr_omega_m = (1 - attempt_progress) * base_omega_m + attempt_progress * target_omega_m;
            p_curr.omega_m = curr_omega_m;
            p_curr.omega_s = 1.0 - curr_omega_m; 
        end
        
        Targets_curr = Targets;
        if do_R0_homotopy
            Targets_curr.R0 = (1 - attempt_progress) * base_R0 + attempt_progress * target_R0_val;
        end
        if do_tech_homotopy
            for i = 1:length(sectors)
                sec = sectors{i};
                p_curr.A0.(sec) = (1 - attempt_progress) * base_tech.A0.(sec) + attempt_progress * target_tech.A0.(sec);
                p_curr.(['g_A_' sec]) = (1 - attempt_progress) * base_tech.g_A.(sec) + attempt_progress * target_tech.g_A.(sec);
            end
        end
        
        Guess.KT = Current_KT; Guess.RT = Current_RT;
        
        fprintf('----------------------------------------------------------\n');
        fprintf('Step %3d: Prog=%6.2f%% (Step: %.1e) | ', ...
                iter_count, attempt_progress*100, step_size);
        print_policy_status(CurrentPolicy);
        fprintf('\n');
        
        try
            [res, best_KT, best_RT, success] = solve_shooting_2sector_inverted(p_curr, CurrentPolicy, Guess, Targets_curr, ShootingSettings);
        catch ME
            fprintf('\n[Error] %s\n', ME.message);
            success = false; best_KT = NaN; best_RT = NaN;
        end
        
        if success
            fprintf('Step converged (KT=%.4f, RT=%.4f)\n', best_KT, best_RT);
            current_progress = attempt_progress;
            Current_KT = best_KT; Current_RT = best_RT;
            FinalRes = res;
            
            History.Progress(end+1) = current_progress;
            History.KT(end+1) = Current_KT; History.RT(end+1) = Current_RT;
            
            step_size = min(step_size * 1.5, max_step_size);
            if current_progress >= 0.99999, break; end
        else
            fprintf('Step failed. Reducing step size (%.2e -> %.2e)\n', step_size, step_size*0.5);
            step_size = step_size * 0.5;
            if step_size < min_step_size
                fprintf('\nMinimum step size reached. Homotopy failed.\n'); 
                break;
            end
        end
    end
    
    if current_progress >= 0.999
        fprintf('\nHomotopy completed successfully.\n');
    else
        fprintf('\nHomotopy incomplete. Final progress: %.2f%%\n', current_progress * 100);
    end
end

function Pol = generate_policy_from_params(Target, alpha, T)
    time_vec = 0:T;
    
    start_t = 0;
    end_t = T;
    if isfield(Target, 'policy_start_t')
        start_t = Target.policy_start_t;
    end
    if isfield(Target, 'policy_end_t')
        end_t = Target.policy_end_t;
    end
    
    mask = (time_vec >= start_t) & (time_vec <= end_t);
    
    Pol.tau_e_path = zeros(1, T + 1);
    if isfield(Target, 'tau_e_0')
        tgt_tau_0 = Target.tau_e_0;
        tgt_g     = 0; 
        if isfield(Target, 'g_tau_e'), tgt_g = Target.g_tau_e; end
        
        curr_tau_0 = tgt_tau_0 * alpha;
        curr_g     = tgt_g * alpha; 
        
        raw_path = curr_tau_0 * ((1 + curr_g) .^ time_vec);
        Pol.tau_e_path(mask) = raw_path(mask);
    end
    
    simple_taxes = {'tau_x', 'tau_int', 'tau_c_m', 'tau_c_s'};
    for i = 1:length(simple_taxes)
        name = simple_taxes{i};
        path_name = [name '_path'];
        
        Pol.(path_name) = zeros(1, T + 1);
        
        if isfield(Target, name)
            tgt_val = Target.(name);
            curr_val = tgt_val * alpha;
            Pol.(path_name)(mask) = curr_val;
        end
    end
end

function print_policy_status(Pol)
    max_tau_e = max(abs(Pol.tau_e_path));
    if max_tau_e > 1e-6
        fprintf('Max Tau_e=%.2f%% ', max_tau_e*100);
    end
    
    max_tau_cs = max(abs(Pol.tau_c_s_path));
    if max_tau_cs > 1e-6
        fprintf('Max Tau_cs=%.2f%% ', max_tau_cs*100);
    end
    
    max_tau_cm = max(abs(Pol.tau_c_m_path));
    if max_tau_cm > 1e-6
        fprintf('Max Tau_cm=%.2f%% ', max_tau_cm*100);
    end
end