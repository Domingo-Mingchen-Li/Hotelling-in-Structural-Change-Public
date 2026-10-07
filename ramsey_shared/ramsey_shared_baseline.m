function [b,history]=ramsey_shared_baseline(base,beta,reference,inherited,w,out,T)
L=inherited.calendar_date;N=T-L;
assert(N>0&&reference.T==T,'RamseyShared:Dates','Invalid calendar horizon');
zero=struct('id','baseline','start',0,'duration',10,'permanent',false, ...
 'tau_c_m',0,'tau_c_s',0,'tau_e',0,'tau_int',0,'tau_x',0,'investment_policy_kind','none');
current=reference;current.T=N;
for f=fieldnames(reference.path).'
 current.path.(f{1})=reference.path.(f{1})(L+1:end);
end
current.Targets.K0=inherited.K;current.Targets.R0=inherited.R;
localbase=base;localbase.targets=current.Targets;localbase.A0=inherited.A;
localbase.scale_h0=current.path.h(1);localbase.scale_E0=current.path.E(1);
localbase.newton_options=w.ramsey.inner_options;
localbase=prepare_config(localbase);
lambda=0;step=w.production_initial_step;history={};
for attempt=1:w.maximum_continuation_attempts
 if lambda>=1,break;end
 next=min(1,lambda+step);cfg=ramsey_shared_production_config(localbase,beta,next);
 cfg.newton_options.display=w.display_production_inner;
 cfg.audit.derivatives=false;cfg.audit.restart=false;

 r=solve_transition(cfg,zero,N,current,[]);
 history{end+1}=struct('attempt',attempt,'lambda',next,'passed',r.passed); 
 if r.passed
  current=r;lambda=next;step=min(w.production_maximum_step,step*w.step_growth);
 else
  step=step/2;
 end
 save(fullfile(out,sprintf('shared_production_checkpoint_T%d.mat',T)),'current','lambda','step','history','beta','inherited');
 assert(step>=w.production_minimum_step,'RamseyShared:Continuation','Continuation below minimum');
end
assert(lambda==1,'RamseyShared:Continuation','Continuation exhausted attempt limit');
cfg=ramsey_shared_production_config(localbase,beta,1);
cfg.newton_options.display=w.display_production_inner;cfg.audit.derivatives=true;cfg.audit.restart=false;
local=solve_transition(cfg,zero,N,current,[]);
assert(local.passed,'RamseyShared:Baseline','Post-switch baseline failed');
state_error=max(abs([local.path.K(1)/inherited.K-1,local.path.R(1)/inherited.R-1]));
assert(state_error<1e-8,'RamseyShared:States','Post-switch initial stocks differ');
V=[local.path.p_m.*local.path.c_m;local.path.p_s.*local.path.c_s;local.path.I];
bar=([local.p.beta_m local.p.beta_s local.p.beta_x]*V)./sum(V,1);
assert(max(abs(bar-beta))<1e-12,'RamseyShared:Channel','Post-switch composition coefficient is not constant');
save(fullfile(out,sprintf('shared_post_switch_baseline_T%d.mat',T)),'local','inherited','state_error');
b=local;b.T=T;b.calendar_terminal=T;b.post_switch_local_horizon=N;
b.history_prefix_is_seed_only=true;b.production_switch_date=L;
for f=fieldnames(local.path).'
 prefix=reference.path.(f{1})(1:L);suffix=local.path.(f{1});
 b.path.(f{1})=[reshape(prefix,1,[]) reshape(suffix,1,[])];
end

end
