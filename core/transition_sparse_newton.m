function sol = transition_sparse_newton(q,x0,options)
opt=struct('max_iterations',50,'residual_tolerance',1e-9,'jacobian_step',1e-5, ...
    'max_log_step',1,'armijo',1e-4,'max_backtracks',25,'minimum_step',1e-12,'display',true);
if nargin>=3&&~isempty(options)
    names=fieldnames(options);
    for k=1:numel(names)
        assert(isfield(opt,names{k}),'v2:Options','Unknown option %s',names{k});
        opt.(names{k})=options.(names{k});
    end
end
assert(opt.max_iterations>=1&&opt.max_iterations==fix(opt.max_iterations) ...
    &&opt.residual_tolerance>0&&opt.jacobian_step>0&&opt.max_log_step>0 ...
    &&opt.armijo>0&&opt.armijo<1&&opt.max_backtracks>=0 ...
    &&opt.max_backtracks==fix(opt.max_backtracks)&&opt.minimum_step>0,'v2:Options','Invalid options');
x=x0(:);sol=struct('success',false,'status','not_started','options',opt, ...
    'x',x,'iterations',0,'residual_evaluations',0,'jacobian_evaluations',0);
[F,d,reason]=evaluate_trial(x,q);sol.residual_evaluations=1;
if ~isempty(reason)
    sol.status=['initial_guess_rejected: ' reason];sol.history=[];return;
end
entry=struct('iteration',0,'maximum_residual',norm(F,inf),'merit',0.5*(F'*F), ...
    'alpha',0,'step_inf',0,'backtracks',0,'feasibility_rejections',0, ...
    'evaluation_rejections',0,'merit_rejections',0,'linear_residual',NaN);
history=entry;
if opt.display,print_entry(entry);end
for k=1:opt.max_iterations
    if norm(F,inf)<=opt.residual_tolerance
        sol.success=true;sol.status='converged';break;
    end
    [J,ji]=transition_sparse_jacobian(x,q,opt.jacobian_step);
    sol.jacobian_evaluations=sol.jacobian_evaluations+1;
    sol.residual_evaluations=sol.residual_evaluations+ji.residual_evaluations;
    direction=-(J\F); 
    if ~isreal(direction)||any(~isfinite(direction))
        sol.status='nonfinite_newton_direction';break;
    end
    linear_error=norm(J*direction+F,inf)/max(norm(F,inf),eps);
    slope=F'*(J*direction);merit=0.5*(F'*F);
    if linear_error>1e-6||~isfinite(slope)||slope>=0
        sol.status='inaccurate_or_non_descent_newton_direction';break;
    end
    alpha=min(1,opt.max_log_step/max(norm(direction,inf),eps));
    accepted=false;infeasible=0;invalid=0;no_decrease=0;
    for bt=0:opt.max_backtracks
        if alpha*norm(direction,inf)<opt.minimum_step,break;end
        trial=x+alpha*direction;
        [Ft,dt,reason]=evaluate_trial(trial,q);
        sol.residual_evaluations=sol.residual_evaluations+1;
        if isempty(reason)
            if 0.5*(Ft'*Ft)<=merit+opt.armijo*alpha*slope
                accepted=true;break;
            end
            no_decrease=no_decrease+1;
        elseif strcmp(reason,'infeasible_or_clamped')
            infeasible=infeasible+1;
        else
            invalid=invalid+1;
        end
        alpha=alpha/2;
    end
    if ~accepted
        sol.status='line_search_failed';
        sol.failed_search=struct('iteration',k,'alpha',alpha,'feasibility_rejections',infeasible, ...
            'evaluation_rejections',invalid,'merit_rejections',no_decrease,'last_reason',reason);
        break;
    end
    entry=struct('iteration',k,'maximum_residual',norm(Ft,inf),'merit',0.5*(Ft'*Ft), ...
        'alpha',alpha,'step_inf',norm(alpha*direction,inf),'backtracks',bt, ...
        'feasibility_rejections',infeasible,'evaluation_rejections',invalid, ...
        'merit_rejections',no_decrease,'linear_residual',linear_error);
    history(end+1)=entry; 
    x=trial;F=Ft;d=dt;sol.iterations=k;
    if opt.display,print_entry(entry);end
end
if norm(F,inf)<=opt.residual_tolerance&&d.all_feasible&&isempty(d.clamp_times)
    sol.success=true;sol.status='converged';
elseif strcmp(sol.status,'not_started')
    sol.status='maximum_iterations_reached';
end
sol.x=x;sol.residual=F;sol.diagnostics=d;sol.history=history;
if opt.display,end
end

function [F,d,reason]=evaluate_trial(x,q)
F=[];d=[];reason='';
try
    [F,d]=transition_global_residual(x,q);
    if ~d.all_feasible||~isempty(d.clamp_times),reason='infeasible_or_clamped';end
catch ME
    numeric_ids={'v2:Coordinates','v2:State','v2:Residual','v2:Return'};
    if any(strcmp(ME.identifier,numeric_ids))
        reason=ME.identifier;
    else
        rethrow(ME); 
    end
end
end

function print_entry(e)

end
