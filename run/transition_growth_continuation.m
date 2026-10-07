function report = transition_growth_continuation(cfg,opt)
report=struct('passed',false,'model_mode',cfg.model,'horizons',cfg.horizons, ...
    'baseline_cases',{{}},'policy_cases',{{}},'warm_starts',{{}},'attempts',{{}}, ...
    'empirical_fit_validated',false,'infinite_horizon_accuracy_validated',false);
try
    raw=load(cfg.warm_start_file,'warm_starts');pool=raw.warm_starts(:).';T=cfg.horizons;
    warm=select_initial_path(pool,cfg.model,'baseline',T);
    fields={'g_A_m','g_A_s','g_A_x','g_A_e'};start=zeros(1,4);target=start;
    original=warm.p;remaining=cfg.p;
    for k=1:4
        start(k)=original.(fields{k});target(k)=remaining.(fields{k});
        original=rmfield(original,fields{k});remaining=rmfield(remaining,fields{k});
    end
    assert(same_run_metadata(original,remaining)&&same_run_metadata(warm.Targets,cfg.targets), ...
        'v2:InitializationScope','This initializer changes ONLY growth. Other calibration/initial stocks differ');
    report.growth_fields=fields;report.start_growth=start;report.target_growth=target;
    zero=struct('id','baseline','permanent',false,'start',0,'duration',11, ...
        'tau_c_m',0,'tau_c_s',0,'tau_int',0,'tau_e',0,'tau_x',0,'investment_policy_kind','none');
    source=cfg;for k=1:4,source.p.(fields{k})=start(k);end
    source.audit.derivatives=false;source.audit.restart=false;
    source.newton_options.display=false;

    solved=solve_transition(source,zero,T,warm,[]);
    report.source_validation=solved;
    assert(solved.passed,'v2:InitializationSource','Cannot validate baseline at SOURCE growth/horizon');
    warm=compact(solved);lambda=0;step=opt.initial_step;
    for attempt=1:opt.maximum_attempts
        if lambda>=1,break;end
        next=min(1,lambda+step);trial_cfg=source;
        g=exp((1-next)*log(start)+next*log(target));
        if next==1,g=target;end 
        for k=1:4,trial_cfg.p.(fields{k})=g(k);end

        trial=solve_transition(trial_cfg,zero,T,warm,[]);
        item=struct('attempt',attempt,'lambda_from',lambda,'lambda_trial',next, ...
            'growth',g,'step',step,'accepted',trial.passed);
        if isfield(trial,'solver'),item.solver_status=trial.solver.status;item.iterations=trial.solver.iterations;end
        if isfield(trial,'seed_diagnostics')
            item.seed_maximum_residual=trial.seed_diagnostics.max_equilibrium_residual;
            item.seed_infeasible_times=trial.seed_diagnostics.infeasible_times;
        end
        if trial.passed
            item.maximum_residual=trial.validation.diagnostics.max_equilibrium_residual;
            report.attempts{end+1}=item;
            warm=compact(trial);lambda=next;
            step=min(opt.maximum_step,step*opt.step_growth);
        else
            item.error_identifier=trial.error_identifier;item.error_message=trial.error_message;
            report.attempts{end+1}=item;report.last_failed_trial=trial;
            assert(any(strcmp(trial.error_identifier,{'v2:Convergence','v2:Tail','v2:Return','v2:State','v2:Coordinates','v2:Residual'})), ...
                'v2:InitializationUnexpected','Unexpected stage error: %s',trial.error_message);
            step=step/2;

            assert(step>=opt.minimum_step,'v2:InitializationStep','Continuation step below minimum');
        end
    end
    report.lambda_reached=lambda;
    assert(lambda>=1,'v2:InitializationAttempts','Target growth not reached within attempt limit');
    final=solve_transition(cfg,zero,T,warm,[]);
    report.baseline_cases={final};
    assert(final.passed,'v2:InitializationFinal','Final current-growth derivative/restart/equilibrium audit failed');
    assert(all(cellfun(@(k) final.p.(k)==cfg.p.(k),fields)), ...
        'v2:InitializationTarget','Final growth factors differ from requested factors');
    report.warm_starts={compact(final)};report.passed=true;
catch ME
    report.error_identifier=ME.identifier;report.error_message=ME.message;

end
end

function s=compact(r)
s=struct('model_mode',r.model_mode,'case_name',r.case_name,'T',r.T,'p',r.p, ...
    'Targets',r.Targets,'policy',r.policy,'investment_policy_kind',r.investment_policy_kind, ...
    'boundary_kind',r.boundary_kind,'path',r.path,'validation',struct('passed',r.validation.passed));
end
