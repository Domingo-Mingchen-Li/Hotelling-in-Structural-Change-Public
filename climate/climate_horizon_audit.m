function a=climate_horizon_audit(r,c,m)
chi=max(abs([c.scenarios.chi]./[r.scenarios.chi]-1));
damage=max(abs([c.scenarios.discounted_damage]./[r.scenarios.discounted_damage]-1));
loss=max(abs([c.scenarios.utility_loss]./[r.scenarios.utility_loss]-1));
v1=[r.inherited.K;r.inherited.R;r.inherited.A;r.inherited.environment_components];
v2=[c.inherited.K;c.inherited.R;c.inherited.A;c.inherited.environment_components];
inherit=max(abs(v2-v1)./max(abs(v1),1e-12));
a=struct('T',r.T,'check_T',c.T,'maximum_chi_relative_difference',chi, ...
 'maximum_damage_relative_difference',damage,'maximum_utility_loss_relative_difference',loss, ...
 'maximum_inherited_state_relative_difference',inherit,'fixed_shared_normalization',true, ...
 'tolerance',m.horizon_relative_tolerance);
a.passed=max([chi damage loss inherit])<=m.horizon_relative_tolerance;

end
