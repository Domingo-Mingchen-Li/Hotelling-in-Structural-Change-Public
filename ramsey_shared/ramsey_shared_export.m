function comparison=ramsey_shared_export(original,flat,meta,out,w)
rows=[];hist=[];comparison=struct('own_history',false,'chi_fixed',true,'normalization_fixed',true, ...
 'common_inherited_states',true,'isolated_causal_channel',false,'global_optimum_certified',false);
for j=1:numel(flat.runs)
 cf=flat.runs{j};T=cf.calendar_terminal;
 idx=find(cellfun(@(a)a.calendar_terminal,original.runs)==T,1);ref=original.runs{idx};
 N=cf.baseline.T;L=flat.common_inherited.calendar_date;
 assert(ref.baseline.T==N&&original.common_inherited.calendar_date==L,'RamseyCF:Dates','Calendar mismatch');
 rp=transition_export_path(ref.best.evaluation.x,ref.best.evaluation.problem);
 cp=transition_export_path(cf.best.evaluation.x,cf.best.evaluation.problem);
 rt=ref.best.evaluation.problem.policy.tau_e_path(:);ct=cf.best.evaluation.problem.policy.tau_e_path(:);
 [rmec,rmuh]=marginals(ref,original.common_normalization.e_reference,meta.climate_kernel_weights);
 [cmec,cmuh]=marginals(cf,flat.common_normalization.e_reference,meta.climate_kernel_weights);
 calendar_date=(L:L+N)';
 tab=table(calendar_date,rt,ct,100*(ct-rt),rmec,cmec,rmuh,cmuh, ...
  100*(rp.R_input(:)./ref.baseline.path.R_input(:)-1), ...
  100*(cp.R_input(:)./cf.baseline.path.R_input(:)-1), ...
  'VariableNames',{'calendar_date','original_optimal_resource_tax','flat_optimal_resource_tax', ...
  'tax_difference_percentage_points','original_marginal_environment_cost','flat_marginal_environment_cost', ...
  'original_mu_times_h','flat_mu_times_h','original_extraction_deviation_percent','flat_extraction_deviation_percent'});
 writetable(tab,fullfile(out,sprintf('optimal_tax_comparison_T%d.csv',T)));
 for k=1:2
  run=ref;name='original';if k==2,run=cf;name='flat_final_input_shares';end
  e=run.best.evaluation;a=e.welfare;b=run.baseline_welfare;
  private=a.private_value-b.private_value;environment=b.discounted_damage-a.discounted_damage;
  row=struct('economy',name,'calendar_terminal',T,'announcement_date',L,'chi',original.scenario.chi, ...
   'initial_resource_tax',e.problem.policy.tau_e_path(1), ...
   'consumption_utility_change',private,'environment_welfare_change',environment,'net_welfare_change',private+environment, ...
   'net_gain_in_fixed_original_anchor_units',e.gain,'projected_gradient',run.best.projected_gradient, ...
   'upper_bound_active',run.upper_bound_active,'own_history',false,'chi_recalibrated',false);
  if isempty(rows),rows=row;else,rows(end+1)=row;end 
 end

end
writetable(struct2table(rows),fullfile(out,'optimal_welfare_comparison.csv'));
for k=1:2
 s=original.common_inherited;name='original';if k==2,s=flat.common_inherited;name='flat_final_input_shares';end
 row=struct('economy',name,'announcement_date',s.calendar_date,'inherited_K',s.K,'inherited_R',s.R, ...
  'inherited_environment_burden',s.environment_burden,'component_1',s.environment_components(1), ...
  'component_2',s.environment_components(2),'component_3',s.environment_components(3),'component_4',s.environment_components(4), ...
  'original_K0',meta.flat_configuration.targets.K0,'original_R0',meta.flat_configuration.targets.R0, ...
  'fixed_e_reference',original.common_normalization.e_reference,'fixed_M_reference',original.common_normalization.M_reference);
 if isempty(hist),hist=row;else,hist(end+1)=row;end 
end
writetable(struct2table(hist),fullfile(out,'inherited_states_comparison.csv'));
comparison.rows=rows;comparison.inherited_states=hist;
comparison.interpretation='Unexpected production switch at the announcement date, conditional on identical inherited economic and climate states';
end
function [mec,muh]=marginals(run,e_ref,a)
e=run.best.evaluation;p=e.problem.p;N=e.problem.T;P=transition_export_path(e.x,e.problem);
mec=nan(N+1,1);mu=P.E.^(p.epsilon-1).*P.P_star.^(-p.epsilon);muh=(mu.*P.h).';
for t=0:N-1,mec(t+1)=p.beta^(-t)*a'*e.welfare.environment_costate_discounted(:,t+2)/e_ref;end
end
