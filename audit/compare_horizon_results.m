function c = compare_horizon_results(a,b,base_a,base_b)
c=struct('T',a.T,'reference_T',b.T,'windows',{{}},'terminal',a.validation.diagnostics.terminal);
for last=unique([10,25,min(a.T,b.T)])
    w=struct('window_end',last,'path',transition_compare_horizons(a,b,last));
    if isfield(a,'effects')
        idx=1:last+1;names=a.effects.variable_names;
        rows=find(~strcmp(names,'E'));delta=a.effects.level_effect_percent(rows,idx)-b.effects.level_effect_percent(rows,idx);
        w.maximum_real_effect_difference_pp=max(abs(delta),[],'all');
        w.real_effect_variable_names=names(rows);w.real_effect_differences_pp=delta;
        w.maximum_output_effect_difference_pp=max(abs(a.effects.output_effect_percent(idx)-b.effects.output_effect_percent(idx)));
        w.maximum_service_share_effect_difference_pp=max(abs(a.effects.service_consumption_share_effect_pp(idx)-b.effects.service_consumption_share_effect_pp(idx)));
        w.maximum_labor_effect_difference_pp=0;
        for sector={'m','s','x','e'}
            key=['l_' sector{1} '_effect_pp'];
            w.maximum_labor_effect_difference_pp=max(w.maximum_labor_effect_difference_pp, ...
                max(abs(a.effects.(key)(idx)-b.effects.(key)(idx))));
        end
    end
    c.windows{end+1}=w;
end
c.early_path=c.windows{find(cellfun(@(w)w.window_end==25,c.windows))}.path;
if isfield(a,'effects')
    ea=a.effects;eb=b.effects;
    c.extraction_t0_effect_difference_pp=ea.extraction_t0_percent-eb.extraction_t0_percent;
    c.cumulative_effect_differences_pp=[ea.extraction_total_0_10_percent-eb.extraction_total_0_10_percent, ...
        ea.extraction_total_0_25_percent-eb.extraction_total_0_25_percent, ...
        ea.extraction_total_0_50_percent-eb.extraction_total_0_50_percent];
    c.cumulative_windows=[10,25,50];
    c.decomposition_component_difference_log_points=a.decomposition.component_log_points-b.decomposition.component_log_points;
    c.total_tilt_difference_log_points=a.decomposition.total_log_points-b.decomposition.total_log_points;
    c.structural_share_difference_pp=NaN;
    if ~a.decomposition.share_denominator_near_zero&&~b.decomposition.share_denominator_near_zero
        c.structural_share_difference_pp=100*(a.decomposition.signed_component_shares(3)-b.decomposition.signed_component_shares(3));
    end
    w=c.windows{find(cellfun(@(w)w.window_end==25,c.windows))};
    c.early_maximum_real_effect_difference_pp=w.maximum_real_effect_difference_pp;
    c.screen_tolerances=struct('early_real_effect_pp',0.05,'early_labor_effect_pp',0.01, ...
        'initial_extraction_pp',0.05,'cumulative_0_25_pp',0.05,'decomposition_log_points',0.05);
    c.screen_checks=[w.maximum_real_effect_difference_pp<=0.05,w.maximum_labor_effect_difference_pp<=0.01, ...
        abs(c.extraction_t0_effect_difference_pp)<=0.05,max(abs(c.cumulative_effect_differences_pp(1:2)))<=0.05, ...
        max(abs(c.decomposition_component_difference_log_points))<=0.05];
    band=1e-6;
    va=[ea.extraction_t0_percent,ea.extraction_total_0_10_percent,ea.extraction_total_0_25_percent,ea.extraction_total_0_50_percent];
    vb=[eb.extraction_t0_percent,eb.extraction_total_0_10_percent,eb.extraction_total_0_25_percent,eb.extraction_total_0_50_percent];
    sa=sign(va);sb=sign(vb);sa(abs(va)<=band)=0;sb(abs(vb)<=band)=0;
    c.extraction_signs=sa;c.reference_extraction_signs=sb;c.extraction_signs_agree=isequal(sa,sb);
    c.baseline_path_comparison=transition_compare_horizons(base_a,base_b,25);
else
    c.screen_tolerances=struct('early_path_percent',0.05,'early_share_pp',0.01);
    c.screen_checks=[100*c.early_path.maximum_over_variables<=0.05, ...
        c.early_path.share_difference_percentage_points<=0.01];
end
c.screen_passed=all(c.screen_checks);
c.infinite_horizon_accuracy_certified=false;
end
