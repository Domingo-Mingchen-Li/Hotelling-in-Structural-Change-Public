function cfg = counterfactual_production_config(base,common_beta,lambda)
assert(isscalar(common_beta)&&isfinite(common_beta)&&common_beta>0&&common_beta<1, ...
 'CF:Shares','Common beta must lie strictly between zero and one');
assert(lambda>=0&&lambda<=1,'CF:Shares','Invalid continuation amplitude');
cfg=base;
for sector={'m','s','x'}
 i=sector{1};a=base.p.(['alpha_' i]);b=base.p.(['beta_' i]);
 newb=(1-lambda)*b+lambda*common_beta;
 cfg.p.(['beta_' i])=newb;
 cfg.p.(['alpha_' i])=(1-newb)*a/(1-b);
end
cfg=prepare_config(cfg);
g=transition_abgp_rates(cfg.p,'paper');
assert(g.consistent&&g.r_ss>0&&g.nonhom_growth_paper<1, ...
 'CF:Tail','Production counterfactual does not support the implemented ABGP candidate');
end
