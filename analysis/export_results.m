function export_results(report,out)
rows={};
for j=1:numel(report.baseline_cases)
    s=report.baseline_cases{j};d=s.validation.diagnostics;
    export_path(s,out);
    writetable(struct2table(s.model_implied_moments),fullfile(out,sprintf('baseline_T%d_moments.csv',s.T)));
    rows(end+1,:)={'baseline',s.T,d.max_equilibrium_residual,NaN,NaN,NaN,NaN,NaN, ...
        NaN,NaN,NaN,NaN,d.terminal.investment_growth_gap,100*max(abs(d.terminal.share_corrections))}; 
end
for k=1:numel(report.policy_cases)
    c=report.policy_cases{k};
    for j=1:numel(c.results)
        s=c.results{j};e=s.effects;d=s.decomposition;
        rows(end+1,:)={c.id,s.T,s.validation.diagnostics.max_equilibrium_residual, ...
            e.extraction_t0_percent,e.policy_window_cumulative_extraction_effect_percent, ...
            e.extraction_total_0_25_percent,e.extraction_total_0_50_percent,d.total_log_points, ...
            d.component_log_points(1),d.component_log_points(2),d.component_log_points(3), ...
            100*d.signed_component_shares(3),s.validation.diagnostics.terminal.investment_growth_gap, ...
            100*max(abs(s.validation.diagnostics.terminal.share_corrections))}; 
        export_path(s,out);
    end
end
t=cell2table(rows,'VariableNames',{'case_id','T','max_residual','initial_extraction_percent', ...
    'policy_window_cumulative_extraction_percent','cumulative_0_25_extraction_percent', ...
    'cumulative_0_50_extraction_percent','total_tilt_log_points','Hotelling_log_points', ...
    'Growth_log_points','Structure_log_points','signed_structure_share_percent', ...
    'terminal_investment_gap_fraction','terminal_share_correction_pp'});
writetable(t,fullfile(out,'policy_summary.csv'));
end

function export_path(s,out)
path=s.path;time=s.time(:);
names={'K','R','r','h','E','I','R_input','c_m','c_s','s_m','s_s','l_m','l_s','l_x','l_e','transfer','investment_user_price'};
values=zeros(numel(time),numel(names));
for n=1:numel(names),values(:,n)=path.(names{n})(:);end
t=array2table([time,values],'VariableNames',[{'time'},names]);
writetable(t,fullfile(out,sprintf('%s_T%d_path.csv',s.case_name,s.T)));
end
