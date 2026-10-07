function cfg = prepare_config(cfg)
assert(any(strcmp(cfg.model,{'paper','legacy'})),'v2:Mode','Unknown model');
assert(isvector(cfg.horizons)&&~isempty(cfg.horizons)&&all(isfinite(cfg.horizons)) ...
    &&all(cfg.horizons==fix(cfg.horizons))&&all(cfg.horizons>=50), ...
    'v2:Horizon','This reusable suite requires integer horizons >=50');
cfg.horizons=unique(cfg.horizons(:).');
p=cfg.p;p.zeta_s=-p.zeta_m; 
assert(abs(p.omega_m+p.omega_s-1)<1e-12&&p.omega_m>0&&p.omega_s>0 ...
    &&p.beta>0&&p.beta<1&&p.delta>=0&&p.delta<1&&p.epsilon>0&&p.xi>0, ...
    'v2:Calibration','Invalid preference/macroeconomic parameters');
assert(cfg.targets.K0>0&&cfg.targets.R0>0&&all(isfinite([cfg.targets.K0,cfg.targets.R0])) ...
    &&numel(cfg.A0)==4&&all(cfg.A0>0)&&all(isfinite(cfg.A0)), ...
    'v2:Calibration','Invalid initial stocks/technology');
for name={'m','s','x','e'}
    k=name{1};a=p.(['alpha_' k]);b=p.(['beta_' k]);g=p.(['g_A_' k]);
    assert(a>0&&b>0&&a+b<1&&g>0&&all(isfinite([a,b,g])),'v2:Technology','Invalid sector %s',k);
    p.sect.(k)=struct('alpha',a,'beta',b);p.A0.(k)=cfg.A0(find(strcmp({'m','s','x','e'},k)));
end
cfg.p=p;
assert(isfinite(cfg.tolerance)&&cfg.tolerance>0&&cfg.tolerance<=1e-9, ...
    'v2:Tolerance','Acceptance tolerance must be positive and <=1e-9');
ids={cfg.policies.id};assert(numel(unique(ids))==numel(ids)&&~any(strcmp(ids,'baseline')), ...
    'v2:Policy','Policy ids must be unique and cannot be baseline');
for i=1:numel(cfg.policies)
    s=cfg.policies(i);
    assert(isvarname(s.id)&&s.start>=0&&s.start==fix(s.start)&&s.duration>=2&&s.duration==fix(s.duration), ...
        'v2:Policy','Invalid policy id/date/duration');
    for name={'tau_c_m','tau_c_s','tau_int','tau_e','tau_x'}
        x=s.(name{1});assert(isscalar(x)&&isreal(x)&&isfinite(x)&&x>-1,'v2:Policy','Invalid tax wedge');
    end
    assert(any(strcmp(s.investment_policy_kind,{'none','purchase_subsidy'})),'v2:Policy','Unknown investment semantics');
    if strcmp(s.investment_policy_kind,'purchase_subsidy')
        assert(strcmp(cfg.model,'paper')&&~s.permanent&&s.start+s.duration-1<min(cfg.horizons)-1, ...
            'v2:Scope','Purchase subsidy must be temporary in paper mode and expire before the tail');
    elseif strcmp(cfg.model,'paper')
        assert(s.tau_x==0,'v2:Scope','Paper tau_x requires purchase_subsidy semantics');
    end
    assert(isscalar(s.permanent)&&any(s.permanent==[false,true]),'v2:Policy','permanent must be a scalar boolean');
    assert(s.start+s.duration-1<=min(cfg.horizons),'v2:Policy','Comparison window exceeds shortest horizon');
    if ~s.permanent
        assert(s.start+s.duration-1<min(cfg.horizons)-1,'v2:Scope','Temporary policy must expire before the terminal tail');
    end
    assert(s.start<min(cfg.horizons),'v2:Policy','Policy begins outside the shortest horizon');
end
end
