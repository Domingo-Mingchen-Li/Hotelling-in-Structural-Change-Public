function [q,x0] = build_transition_problem(cfg,spec,T,warm)
p=cfg.p;g=transition_abgp_rates(p,'paper');nonhom=g.nonhom_growth_paper;
if strcmp(cfg.model,'legacy'),nonhom=g.nonhom_growth_legacy;end
assert(g.consistent&&g.r_ss>0&&nonhom<1,'v2:Tail','Parameters do not support the implemented ABGP candidate');
q=struct('mode',cfg.model,'p',p,'T',T,'time',0:T,'n',5*(T+1), ...
    'targets',cfg.targets,'growth',g,'nonhom_growth',nonhom, ...
    'boundary_kind','full_static_abgp_tail','investment_policy_kind',spec.investment_policy_kind, ...
    'case_name',spec.id,'policy_specification',spec);
rates=[p.g_A_m;p.g_A_s;p.g_A_x;p.g_A_e];q.technology=cfg.A0(:).*rates.^q.time;
mask=q.time>=spec.start&q.time<=spec.start+spec.duration-1;
if spec.permanent,mask=q.time>=spec.start;end
for name={'tau_c_m','tau_c_s','tau_int','tau_e','tau_x'},q.policy.([name{1} '_path'])=spec.(name{1})*double(mask);end
q.reference=[cfg.targets.K0*g.g_K.^q.time;cfg.targets.R0*g.g_e.^q.time; ...
    g.r_ss*ones(1,T+1);cfg.scale_h0*g.g_h.^q.time;cfg.scale_E0*g.g_E.^q.time];
assert(all(isfinite(q.reference(:)))&&all(q.reference(:)>0),'v2:Scale','Invalid normalized reference');
q.variable_names={'K','R','r','h','E'};
q.equation_names={'capital','capital_accumulation','resource_accumulation','Euler','Hotelling', ...
    'initial_K','initial_R','terminal_r','terminal_capital','terminal_resource_tail'};
X=zeros(5,warm.T+1);for k=1:5,X(k,:)=warm.path.(q.variable_names{k})(:).';end
X(1,:)=X(1,:)*(cfg.targets.K0/warm.Targets.K0);
X(2,:)=X(2,:)*(cfg.targets.R0/warm.Targets.R0);
if T<=warm.T,X=X(:,1:T+1);else,X=[X,X(:,end).*[g.g_K;g.g_e;1;g.g_h;g.g_E].^(1:T-warm.T)];end
q.seed=X;q.seed_strategy='same-model validated path, calendar truncation/extension and initial-stock scaling';
q.groups=cell(1,10);for k=1:5,q.groups{k}=k:5:5*T;end
for k=6:10,q.groups{k}=5*T+k-5;end
q.row_names=cell(q.n,1);
for t=0:T-1,for k=1:5,q.row_names{5*t+k}=sprintf('%s[t=%d]',q.equation_names{k},t);end,end
for k=6:10,q.row_names{5*T+k-5}=q.equation_names{k};end
[q.pattern,q.colors]=transition_jacobian_structure(T);x0=transition_pack_path(X,q);
end
