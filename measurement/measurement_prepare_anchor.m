function measurement_prepare_anchor(m,c,case_root)
seed_file=m.seed_file;
if isfield(m,'calibration_seed_file')&&~isempty(m.calibration_seed_file)
 seed_file=m.calibration_seed_file;
end
assert(isfile(seed_file),'measurement:Seed','Seed file missing: %s',seed_file);
raw=load(seed_file,'warm_starts');w=select_initial_path(raw.warm_starts,'paper','baseline',c.T);
zero=struct('id','baseline','permanent',false,'start',0,'duration',10, ...
 'tau_c_m',0,'tau_c_s',0,'tau_int',0,'tau_e',0,'tau_x',0,'investment_policy_kind','none');
cfg=c.base;cfg.audit.derivatives=false;cfg.audit.restart=false;
cfg.targets=w.Targets;cfg.p.zeta_m=w.p.zeta_m;cfg.p.zeta_s=-cfg.p.zeta_m;
cfg.policy_ids={};cfg.policies=cfg.policies([]);cfg=prepare_config(cfg);
attempts={};r=solve_transition(cfg,zero,c.T,w,[]);
attempts{end+1}=struct('lambda',1,'passed',r.passed);
if ~r.passed
 target=cfg;origin=cfg;
 for sector={'m','s','x','e'}
  for kind={'alpha','beta'}
   key=[kind{1} '_' sector{1}];origin.p.(key)=w.p.(key);
  end
 end
 origin=prepare_config(origin);r=solve_transition(origin,zero,c.T,w,[]);
 assert(r.passed,'measurement:Anchor','Cannot revalidate seed baseline at its production shares');
 lambda=0;step=m.production_continuation_step;
 assert(step>0&&step<=1&&m.production_continuation_minimum_step>0&& ...
 m.production_continuation_minimum_step<=step,'measurement:Options','Invalid continuation steps');
 while lambda<1
  next=min(1,lambda+step);stage=target;
  for sector={'m','s','x','e'}
   for kind={'alpha','beta'}
    key=[kind{1} '_' sector{1}];stage.p.(key)=origin.p.(key)+next*(target.p.(key)-origin.p.(key));
   end
  end
  stage=prepare_config(stage);trial=solve_transition(stage,zero,c.T,r,[]);
  attempts{end+1}=struct('lambda',next,'passed',trial.passed); 
  if trial.passed
   lambda=next;r=trial;step=min(2*step,1-lambda);
  else
   step=step/2;
   assert(step>=m.production_continuation_minimum_step,'measurement:Anchor', ...
    'Production continuation failed; use a better validated baseline seed');
  end
 end
end
folder=fullfile(case_root,'anchor');if ~isfolder(folder),mkdir(folder);end
warm_starts={r};save(fullfile(folder,'initial_paths.mat'),'warm_starts');
save(fullfile(folder,'anchor_audit.mat'),'attempts');
end
