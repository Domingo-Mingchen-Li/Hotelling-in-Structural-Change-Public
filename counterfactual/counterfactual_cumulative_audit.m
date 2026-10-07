function rows = counterfactual_cumulative_audit(model,m)
rows=repmat(struct('case_id','','H',NaN,'T',NaN,'check_T',NaN, ...
 'effect_percent',NaN,'check_effect_percent',NaN,'difference_pp',NaN,'passed',false),1,0);
for k=1:numel(model.policy_cases)
 a=model.policy_cases{k}.results{1};b=model.policy_cases{k}.results{end};
 ba=model.baseline_cases{1};bb=model.baseline_cases{end};
 for H=m.cumulative_periods
  va=100*(sum(a.path.R_input(1:H))/sum(ba.path.R_input(1:H))-1);
  vb=100*(sum(b.path.R_input(1:H))/sum(bb.path.R_input(1:H))-1);
  difference=va-vb;ok=abs(difference)<=m.horizon_cumulative_tolerance_pp;
  rows(end+1)=struct('case_id',a.case_name,'H',H,'T',a.T,'check_T',b.T, ...
   'effect_percent',va,'check_effect_percent',vb,'difference_pp',difference,'passed',ok); 

 end
end
end
