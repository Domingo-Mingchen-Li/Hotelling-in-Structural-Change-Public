function [anchor,history] = counterfactual_prepare_anchor(base,beta,reference,m,out)
zero=struct('id','baseline','permanent',false,'start',0,'duration',10, ...
 'tau_c_m',0,'tau_c_s',0,'tau_int',0,'tau_e',0,'tau_x',0,'investment_policy_kind','none');
current=reference;lambda=0;step=m.production_initial_step;history={};
trial_cfg=base;trial_cfg.audit.derivatives=false;trial_cfg.audit.restart=false;
for attempt=1:m.maximum_continuation_attempts
 if lambda>=1,break;end
 next=min(1,lambda+step);cfg=counterfactual_production_config(trial_cfg,beta,next);
 g=transition_abgp_rates(cfg.p,'paper');
 drawnow limitrate;
 clock=tic;r=solve_transition(cfg,zero,m.T,current,[]);
 history{end+1}=struct('attempt',attempt,'lambda',next,'passed',r.passed,'elapsed',toc(clock)); 
 if r.passed
  lambda=next;current=r;step=min(m.production_maximum_step,step*m.step_growth);

 else
  step=step/2;
 end
 save(fullfile(out,'production_checkpoint.mat'),'current','lambda','step','history','m','beta');
 assert(step>=m.production_minimum_step,'CF:Continuation', ...
  'Production step below minimum at lambda %.8f. Inspect production_checkpoint.mat and console.txt',lambda);
end
assert(lambda==1,'CF:Continuation','Production continuation exceeded attempt limit');
cfg=counterfactual_production_config(base,beta,1);

anchor=solve_transition(cfg,zero,m.T,current,[]);
assert(anchor.passed,'CF:Anchor','Final flat-share baseline audit failed');
warm_starts={anchor};save(fullfile(out,'counterfactual_initial_paths.mat'),'warm_starts');
end
