function r = solve_transition(cfg,spec,T,warm,baseline)
r=struct('passed',false,'case_name',spec.id,'model_mode',cfg.model,'T',T);
try
    [q,x0]=build_transition_problem(cfg,spec,T,warm);r.problem=q;
    [~,r.seed_diagnostics]=transition_global_residual(x0,q);
    sol=transition_sparse_newton(q,x0,cfg.newton_options);
    r.continuation=struct('used',false,'direct_attempt_status',sol.status,'stages',{{}});
    if ~sol.success&&~isempty(baseline)
        r.continuation.used=true;
        X=zeros(5,T+1);for k=1:5,X(k,:)=baseline.path.(q.variable_names{k})(:).';end
        x=transition_pack_path(X,q);
        for lambda=[0.25,0.5,0.75,1]
            qi=q;for name=fieldnames(q.policy).',qi.policy.(name{1})=lambda*q.policy.(name{1});end
            si=transition_sparse_newton(qi,x,cfg.newton_options);vi=transition_validate_solution(si.x,qi,cfg.tolerance);
            r.continuation.stages{end+1}=struct('lambda',lambda,'solver',si,'validation',vi);
            assert(si.success&&vi.passed,'v2:Continuation','Policy continuation failed at %.2f',lambda);
            x=si.x;sol=si;
        end
    end
    r.solver=sol;r.x=sol.x;r.reference=q.reference;
    r.validation=transition_validate_solution(sol.x,q,cfg.tolerance);
    assert(sol.success&&r.validation.passed,'v2:Convergence','Equilibrium failed: %s',sol.status);
    r.path=transition_export_path(sol.x,q);
    if isempty(baseline),r.model_implied_moments=transition_baseline_moments(r.path);end
    r.policy_equation_audit=audit_policy_equations(r.path,q);
    r.p=q.p;r.Targets=q.targets;r.policy=q.policy;r.boundary_kind=q.boundary_kind;
    r.investment_policy_kind=q.investment_policy_kind;r.policy_specification=spec;r.time=q.time;
    if ~isempty(baseline)
        dates=[spec.start,spec.start+spec.duration-1];
        assert(dates(2)>dates(1)&&dates(2)<=T,'v2:Decomposition','Decomposition needs >=2 active nodes within horizon');
        r.decomposition=transition_resource_decomposition(r.path,q.policy,baseline.path,baseline.policy,q.p,dates);
        r.effects=transition_policy_effects(r.path,baseline.path,r.decomposition);
        idx=dates(1)+1:dates(2)+1;
        r.effects.policy_window_dates=dates;
        r.effects.policy_window_cumulative_extraction_effect_percent=100*(sum(r.path.R_input(idx))/sum(baseline.path.R_input(idx))-1);
    end
    if cfg.audit.derivatives,r.derivative_audit=audit_sparse_derivatives(sol.x,q);end
    if cfg.audit.restart
        sp=transition_sparse_newton(q,sol.x+0.001*sin((1:q.n)'*0.73),struct('display',false));
        vp=transition_validate_solution(sp.x,q,cfg.tolerance);err=max(abs(expm1(sp.x-sol.x)));
        r.restart_audit=struct('passed',sp.success&&vp.passed&&err<1e-6,'path_difference',err,'validation',vp);
        assert(r.restart_audit.passed,'v2:Restart','Deterministic local restart failed');
    end
    same=warm.T==T&&strcmp(warm.model_mode,cfg.model)&&strcmp(warm.case_name,spec.id) ...
        &&same_run_metadata(warm.p,cfg.p)&&same_run_metadata(warm.Targets,cfg.targets) ...
        &&same_run_metadata(warm.policy,q.policy)&&strcmp(warm.investment_policy_kind,q.investment_policy_kind);
    r.regression=struct('applicable',same,'maximum_relative_difference',NaN,'metadata_match_ulp_tolerance',32);
    if same
        err=0;for name={'K','R','r','h','E','I','R_input','c_m','c_s'}
            err=max(err,max(abs(r.path.(name{1})./warm.path.(name{1})-1)));
        end
        r.regression.maximum_relative_difference=err;
        assert(err<1e-6,'v2:Regression','Unchanged validated equilibrium changed');
    end
    if strcmp(cfg.model,'paper')&&spec.permanent&&spec.tau_c_m==spec.tau_c_s ...
            &&spec.start==0&&spec.tau_e==0&&spec.tau_int==0&&spec.tau_x==0&&~isempty(baseline)
        err=max(abs(r.path.E./baseline.path.E/(1+spec.tau_c_m)-1));
        for name={'K','R','r','h','I','R_input','c_m','c_s','w','p_m','p_s'}
            err=max(err,max(abs(r.path.(name{1})./baseline.path.(name{1})-1)));
        end
        r.neutrality_error=err;assert(err<1e-6,'v2:Neutrality','Uniform-tax neutrality failed');
    end
    r.passed=true;

catch ME
    r.error_identifier=ME.identifier;r.error_message=ME.message;

end
end
