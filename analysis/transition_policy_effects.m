function e = transition_policy_effects(path,base,decomposition)
names={'K','R','r','h','E','I','R_input','c_m','c_s','w','p_m','p_s'};
e=struct('variable_names',{names},'level_effect_percent',zeros(numel(names),numel(path.K)));
for k=1:numel(names)
    a=path.(names{k})(:).';b=base.(names{k})(:).';
    e.level_effect_percent(k,:)=100*(a./b-1);
end
a=decomposition.policy_objects;b=decomposition.baseline_objects;
e.output_effect_percent=100*(a.Y./b.Y-1);
e.resource_intensity_effect_percent=100*(a.resource_intensity./b.resource_intensity-1);
e.beta_bar_effect_percent=100*(a.beta_bar./b.beta_bar-1);
e.extraction_t0_percent=100*(a.e(1)/b.e(1)-1);
for last=[10,25,50]
    e.(sprintf('extraction_total_0_%d_percent',last))=100*(sum(a.e(1:last+1))/sum(b.e(1:last+1))-1);
end
e.service_consumption_share_effect_pp=100*(path.s_s(:).'-base.s_s(:).');
for name={'m','s','x','e'}
    key=['l_' name{1}];e.([key '_effect_pp'])=100*(path.(key)(:).'-base.(key)(:).');
end
e.permanent_unextracted_stock_validated=false;
end
