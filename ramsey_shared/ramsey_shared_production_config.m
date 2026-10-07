function cfg=ramsey_shared_production_config(base,beta,lambda)
assert(isscalar(beta)&&isfinite(beta)&&beta>0&&beta<1&&lambda>=0&&lambda<=1, ...
 'RamseyCF:Shares','Invalid common beta or continuation amplitude');
cfg=base;
for sector={'m','s','x'}
 i=sector{1};a=base.p.(['alpha_' i]);b=base.p.(['beta_' i]);newb=(1-lambda)*b+lambda*beta;
 cfg.p.(['beta_' i])=newb;cfg.p.(['alpha_' i])=(1-newb)*a/(1-b);
end
cfg=prepare_config(cfg);g=transition_abgp_rates(cfg.p,'paper');
assert(g.consistent&&g.r_ss>0&&g.nonhom_growth_paper<1&&g.g_e<1, ...
 'RamseyCF:Tail','Flattened technology does not support the implemented tail');
assert(cfg.p.beta*(g.g_E/g.g_P)^cfg.p.epsilon<1&& ...
 cfg.p.beta*g.g_pm^cfg.p.zeta_m*g.g_ps^cfg.p.zeta_s<1, ...
 'RamseyCF:Utility','Discounted utility tail diverges');
end
