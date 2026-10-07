function m=climate_prepare_config(m)
assert(strcmp(m.base.model,'paper'),'Climate:Mode','Use paper mode');
assert(isscalar(m.T)&&m.T==fix(m.T)&&m.T>=50&& ...
 isscalar(m.check_T)&&m.check_T==fix(m.check_T)&&m.check_T>m.T,'Climate:Horizon','Invalid T/check_T');
assert(isscalar(m.L)&&m.L==fix(m.L)&&m.L>=1&&m.L<m.T,'Climate:History','Require 1 <= L < T');
assert(isscalar(m.reference_date)&&m.reference_date==fix(m.reference_date)&&m.reference_date>=m.L&&m.reference_date<m.T, ...
 'Climate:Reference','reference_date must be between L and T-1');
m.a=m.a(:);m.timescales_years=m.timescales_years(:);m.initial_components=m.initial_components(:);
assert(numel(m.a)==4&&all(isfinite(m.a))&&all(m.a>=0)&&abs(sum(m.a)-1)<1e-12, ...
 'Climate:Kernel','Four nonnegative weights must sum to one');
assert(numel(m.timescales_years)==4&&all(m.timescales_years>0)&& ...
 all(isfinite(m.timescales_years)|isinf(m.timescales_years)),'Climate:Kernel','Invalid timescales');
assert(numel(m.initial_components)==4&&all(isfinite(m.initial_components))&&all(m.initial_components>=0), ...
 'Climate:Initial','Invalid initial environmental components');
assert(all(ismember(m.etas,[1 2]))&&numel(unique(m.etas))==numel(m.etas),'Climate:Damage','This exact tail supports eta=1 or 2');
assert(all(isfinite(m.loss_targets))&&all(m.loss_targets>0)&all(m.loss_targets<0.2)&& ...
 numel(unique(m.loss_targets))==numel(m.loss_targets),'Climate:Loss','Invalid loss targets');
assert(ismember(m.benchmark_eta,m.etas)&&ismember(m.benchmark_loss,m.loss_targets),'Climate:Benchmark','Benchmark missing');
assert(any(strcmp(m.loss_definition,{'proportional_basket','fixed_price_expenditure'})),'Climate:CE','Unknown definition');
assert(m.maximum_utility_tail_years>=2*m.utility_tail_block_years&& ...
 m.maximum_utility_tail_years==fix(m.maximum_utility_tail_years)&& ...
 m.utility_tail_block_years>=10&&m.utility_tail_block_years==fix(m.utility_tail_block_years), ...
 'Climate:Options','Invalid tail iteration limits');
for f={'utility_tail_relative_tolerance','horizon_relative_tolerance','identity_tolerance','baseline_tail_gap_tolerance'}
 assert(isscalar(m.(f{1}))&&isfinite(m.(f{1}))&&m.(f{1})>0,'Climate:Options','Invalid tolerance');
end
m.rho=exp(-1./m.timescales_years);
m.base.policies=m.base.policies([]);m.base.policy_ids={};
m.base.horizons=[m.T m.check_T];m.base=prepare_config(m.base);
g=transition_abgp_rates(m.base.p,'paper');
gu=(g.g_E/g.g_P)^m.base.p.epsilon;
gpi=g.g_pm^m.base.p.zeta_m*g.g_ps^m.base.p.zeta_s;
assert(m.base.p.beta*gu<1&&m.base.p.beta*gpi<1,'Climate:Utility','Discounted utility tail does not converge');
end
