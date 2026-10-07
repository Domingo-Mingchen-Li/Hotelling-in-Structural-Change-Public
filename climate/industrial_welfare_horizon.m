function h=industrial_welfare_horizon(main,check,anchor,w)
assert(numel(main)==numel(check),'Welfare:Horizon','Missing comparison rows');
h=[];
for k=1:numel(main)
 a=main(k);b=check(k);
 assert(strcmp(a.policy,b.policy)&&a.eta==b.eta&&a.loss_target==b.loss_target,'Welfare:Horizon','Row mismatch');
 errs=abs([a.consumption_utility_change-b.consumption_utility_change, ...
  a.environment_welfare_change-b.environment_welfare_change,a.net_welfare_change-b.net_welfare_change])/anchor;
 same=[a.consumption_sign a.environment_sign a.net_sign]==[b.consumption_sign b.environment_sign b.net_sign];
 r=struct('policy',a.policy,'eta',a.eta,'loss_target',a.loss_target, ...
  'consumption_difference_in_anchor_units',errs(1),'environment_difference_in_anchor_units',errs(2), ...
  'net_difference_in_anchor_units',errs(3),'consumption_sign_stable',same(1), ...
  'environment_sign_stable',same(2),'net_sign_stable',same(3), ...
  'passed',all(errs<w.horizon_anchor_tolerance)&&all(same));
 if isempty(h),h=r;else,h(end+1)=r;end 
end

end
