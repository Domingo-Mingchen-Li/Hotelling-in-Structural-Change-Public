function report = counterfactual_run_cases(cfg,pool,m,out,label)
report=struct('passed',false,'model_mode',cfg.model,'horizons',cfg.horizons, ...
 'baseline_cases',{{}},'policy_cases',{{}},'configuration',cfg, ...
 'empirical_fit_validated',false,'infinite_horizon_accuracy_validated',false);
zero=struct('id','baseline','permanent',false,'start',0,'duration',10, ...
 'tau_c_m',0,'tau_c_s',0,'tau_int',0,'tau_e',0,'tau_x',0,'investment_policy_kind','none');
for j=1:numel(cfg.horizons)
 T=cfg.horizons(j);drawnow limitrate;
 warm=select_initial_path(pool,cfg.model,'baseline',T);
 r=solve_transition(cfg,zero,T,warm,[]);report.baseline_cases{j}=r;
 save(fullfile(out,[label '_checkpoint.mat']),'report');
 assert(r.passed,'CF:Baseline','%s baseline T%d failed; checkpoint saved',label,T);
 pool{end+1}=r; 
end
for k=1:numel(cfg.policies)
 spec=cfg.policies(k);record=struct('id',spec.id,'specification',spec,'results',{{}},'passed',false);
 for j=1:numel(cfg.horizons)
  T=cfg.horizons(j);drawnow limitrate;
  warm=select_initial_path(pool,cfg.model,spec.id,T);
  [r,stages]=policy_solve(cfg,spec,T,warm,report.baseline_cases{j},m);
  r.adaptive_policy_continuation=stages;record.results{j}=r;
  report.policy_cases{k}=record;save(fullfile(out,[label '_checkpoint.mat']),'report');
  assert(r.passed,'CF:Policy','%s %s T%d failed; inspect saved checkpoint',label,spec.id,T);
  pool{end+1}=r; 
 end
 record.passed=true;report.policy_cases{k}=record;
end
report.passed=true;
if numel(cfg.horizons)>1,report=audit_horizons(report);end
warm_starts=pool;save(fullfile(out,[label '_initial_paths.mat']),'warm_starts');
save(fullfile(out,[label '_checkpoint.mat']),'report');
end

function [r,stages]=policy_solve(cfg,spec,T,warm,baseline,m)
clock=tic;r=solve_transition(cfg,spec,T,warm,baseline);stages={};

if r.passed,return;end

current=baseline;lambda=0;step=m.policy_initial_step;
for attempt=1:m.maximum_continuation_attempts
 next=min(1,lambda+step);stage=spec;
 for key={'tau_c_m','tau_c_s','tau_int','tau_e','tau_x'},stage.(key{1})=next*spec.(key{1});end
 drawnow limitrate;
 trial=solve_transition(cfg,stage,T,current,baseline);
 stages{end+1}=struct('lambda',next,'passed',trial.passed); 
 if trial.passed
  lambda=next;current=trial;step=min(1,step*m.step_growth);
  if lambda==1,r=current;return;end
 else
  step=step/2;
  if step<m.policy_minimum_step,r=trial;return;end
 end
end
r=current;r.passed=false;r.error_message='Adaptive policy continuation attempt limit';
end
