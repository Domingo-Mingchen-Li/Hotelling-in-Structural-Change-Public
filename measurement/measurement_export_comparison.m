function comparison = measurement_export_comparison(ref,alt,m,id)
T=m.T;Hs=m.cumulative_periods;
assert(all(isfinite(Hs))&&all(Hs==fix(Hs))&&all(Hs>0)&&all(Hs<=T), ...
 'measurement:Window','Cumulative periods must be positive integers <=T');
assert(numel(unique(Hs))==numel(Hs),'measurement:Window','Duplicate cumulative window');
rows=repmat(struct('measurement_case','','policy','','metric','', ...
 'reference',NaN,'alternative',NaN,'difference',NaN),1,0);
gaps=repmat(struct('policy','','periods',NaN,'relative_identity_error',NaN),1,0);
out=alt.output_dir;
for k=1:numel(alt.policy_cases)
 a=alt.policy_cases{k}.results{1};pid=alt.policy_cases{k}.id;
 j=find(cellfun(@(r)strcmp(r.id,pid),ref.policy_cases));
 assert(numel(j)==1,'measurement:Compare','Missing matching reference policy');
 b=ref.policy_cases{j}.results{1};ab=alt.baseline_cases{1};bb=ref.baseline_cases{1};
 assert(isequaln(a.policy_specification,b.policy_specification), ...
  'measurement:Compare','Reference and alternative policies differ');
 va=metrics(a,ab,Hs,T);vb=metrics(b,bb,Hs,T);names=fieldnames(va);
 for n=1:numel(names)
  key=names{n};r=struct('measurement_case',id,'policy',pid,'metric',key, ...
   'reference',vb.(key),'alternative',va.(key),'difference',va.(key)-vb.(key));
  rows(end+1)=r; 
 end
 for H=Hs
  da=sum(a.path.R_input(1:H))-sum(ab.path.R_input(1:H));
  db=sum(b.path.R_input(1:H))-sum(bb.path.R_input(1:H));
  sa=ab.path.R(H+1)-a.path.R(H+1);sb=bb.path.R(H+1)-b.path.R(H+1);
  err=max(abs([da-sa,db-sb])./max(1,abs([sa,sb])));
  assert(err<1e-5,'measurement:Cumulative','Resource stock/flow identity failed');
  gaps(end+1)=struct('policy',pid,'periods',H,'relative_identity_error',err); 
 end
end
if ~isempty(rows),writetable(struct2table(rows),fullfile(out,'measurement_comparison.csv'));end
pa=alt.configuration.p;pb=ref.configuration.p;
parameters=repmat(struct('parameter','','reference',NaN,'alternative',NaN),1,0);
keys={'alpha_m','beta_m','alpha_s','beta_s','alpha_x','beta_x','alpha_e','beta_e','zeta_m'};
for j=1:numel(keys)
 key=keys{j};parameters(end+1)=struct('parameter',key,'reference',pb.(key),'alternative',pa.(key)); 
end
for key={'K0','R0'}
 parameters(end+1)=struct('parameter',key{1},'reference',ref.configuration.targets.(key{1}), ...
  'alternative',alt.configuration.targets.(key{1})); 
end
writetable(struct2table(parameters),fullfile(out,'measurement_parameters_comparison.csv'));
comparison=struct('rows',rows,'parameters',parameters,'stock_flow_checks',gaps, ...
 'cumulative_periods',Hs,'cumulative_convention','H flows, t=0,...,H-1', ...
 'decomposition_units','log points; exact equilibrium accounting, not causal effects');
save(fullfile(out,'measurement_comparison.mat'),'comparison');
end

function v=metrics(r,base,Hs,T)
d=r.decomposition;obj=d.policy_objects;bobj=d.baseline_objects;
active=r.policy_specification.start+(1:r.policy_specification.duration);
v.extraction_initial_percent=100*(obj.e(1)/bobj.e(1)-1);
v.intensity_initial_percent=100*(obj.resource_intensity(1)/bobj.resource_intensity(1)-1);
v.K_initial_percent=100*(r.path.K(1)/base.path.K(1)-1);
v.K_policy_end_percent=100*(r.path.K(active(end))/base.path.K(active(end))-1);
v.intensity_policy_window_mean_percent=mean(100*(obj.resource_intensity(active)./bobj.resource_intensity(active)-1));
v.intensity_policy_window_max_percent=max(100*(obj.resource_intensity(active)./bobj.resource_intensity(active)-1));
v.extraction_policy_window_cumulative_percent=100*(sum(obj.e(active))/sum(bobj.e(active))-1);
v.total_log_points=d.total_log_points;
v.hotelling_log_points=d.component_log_points(1);
v.growth_log_points=d.component_log_points(2);
v.structural_log_points=d.component_log_points(3);
v.structural_signed_share=d.signed_component_shares(3); 
for H=Hs
 v.(sprintf('cumulative_%d_periods_percent',H))=100*(sum(obj.e(1:H))/sum(bobj.e(1:H))-1);
end
cum=100*(cumsum(obj.e(1:T))./cumsum(bobj.e(1:T))-1);
v.cumulative_min_percent=min(cum);v.cumulative_max_percent=max(cum);
v.cumulative_positive_windows=sum(cum>1e-6);
v.cumulative_negative_windows=sum(cum< -1e-6);
v.cumulative_neutral_windows=sum(abs(cum)<=1e-6);
end

