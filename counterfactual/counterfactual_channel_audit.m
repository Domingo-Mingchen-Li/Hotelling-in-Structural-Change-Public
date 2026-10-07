function [rows,passed] = counterfactual_channel_audit(model,beta,m)
rows=repmat(struct('case_id','','T',NaN,'beta_bar_relative_error',NaN, ...
 'structural_log_points',NaN,'maximum_date_structural_log_points',NaN,'passed',false),1,0);
for j=1:numel(model.horizons)
 b=model.baseline_cases{j};V=[b.path.p_m(:).*b.path.c_m(:),b.path.p_s(:).*b.path.c_s(:),b.path.I(:)];
 beta_bar=V*[b.p.beta_m;b.p.beta_s;b.p.beta_x]./sum(V,2);
 err=max(abs(beta_bar/beta-1));
 rows(end+1)=struct('case_id','baseline','T',b.T,'beta_bar_relative_error',err, ...
  'structural_log_points',0,'maximum_date_structural_log_points',0,'passed',err<=m.beta_bar_relative_tolerance); 
 for k=1:numel(model.policy_cases)
  r=model.policy_cases{k}.results{j};d=r.decomposition;
  err=max(abs(d.policy_objects.beta_bar/beta-1));z=d.component_log_points(3);
  whole=100*max(abs(d.initial_to_each_date_components(3,:)));
  ok=err<=m.beta_bar_relative_tolerance&&abs(z)<=m.structural_zero_tolerance_log_points&& ...
   whole<=m.structural_zero_tolerance_log_points;
  rows(end+1)=struct('case_id',r.case_name,'T',r.T,'beta_bar_relative_error',err, ...
   'structural_log_points',z,'maximum_date_structural_log_points',whole,'passed',ok); 

 end
end
passed=all([rows.passed]);
end
