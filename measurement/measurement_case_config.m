function [cfg,c,spec,signature,case_root] = measurement_case_config(m,id)
assert(isvarname(id),'measurement:Case','Invalid case id');
ids={m.cases.id};assert(numel(unique(ids))==numel(ids),'measurement:Case','Duplicate case ids');
k=find(strcmp(ids,id));assert(numel(k)==1,'measurement:Case','Unknown measurement case %s',id);
spec=m.cases(k);A=spec.shares;
assert(isequal(size(A),[4,2])&&isreal(A)&&all(isfinite(A(:)))&& ...
 all(A(:)>0)&&all(sum(A,2)<1),'measurement:Shares','Expected feasible 4x2 alpha/beta matrix');
assert(isscalar(m.T)&&m.T>=50&&m.T==fix(m.T),'measurement:Horizon','Invalid T');
cfg=m.base;cfg.model='paper';cfg.horizons=m.T;
for j=1:4
 names={'m','s','x','e'};cfg.p.(['alpha_' names{j}])=A(j,1);cfg.p.(['beta_' names{j}])=A(j,2);
end
cfg.policy_ids=m.policy_ids;
assert(all(ismember(m.policy_ids,{cfg.policies.id})),'measurement:Policy','Policy id missing from solver_config active list');
cfg.policies=cfg.policies(ismember({cfg.policies.id},m.policy_ids));
cfg.make_plots=m.make_plots;cfg.audit.derivatives=m.audit_derivatives;
cfg.audit.restart=m.audit_restart;cfg.audit.horizons=false;
cfg=prepare_config(cfg);
case_root=fullfile(m.output_root,id);
c=m.calibration;c.base=cfg;c.T=m.T;
c.base.warm_start_file=fullfile(case_root,'anchor','initial_paths.mat');
if ~isempty(spec.targets),c.targets=spec.targets(:);end
c.output_root=fullfile(case_root,'calibration');
assert(~isempty(c.check_T),'measurement:Calibration','Keep check_T enabled for calibration acceptance');
p=cfg.p;
for name={'zeta_m','zeta_s','sect','A0'}
 if isfield(p,name{1}),p=rmfield(p,name{1});end
end
settings=rmfield(c,{'base','output_root'});
signature=struct('id',id,'shares',A,'fixed_p',p,'A0',cfg.A0, ...
 'fit_settings',settings,'model','paper');
end
