function exports = counterfactual_export(report,out,m)
folder=fullfile(out,'exports');if ~isfolder(folder),mkdir(folder);end
models={report.reference,report.flat};labels={'Original','Equal resource shares'};
keys={'alpha_m','beta_m','alpha_s','beta_s','alpha_x','beta_x','alpha_e','beta_e', ...
 'zeta_m','K0','R0','g_A_m','g_A_s','g_A_x','g_A_e'};values=zeros(numel(keys),2);
for j=1:2
 cfg=models{j}.configuration;
 for k=1:numel(keys)
  if ismember(keys{k},{'K0','R0'}),values(k,j)=cfg.targets.(keys{k});else,values(k,j)=cfg.p.(keys{k});end
 end
end
parameters=table(string(keys(:)),values(:,1),values(:,2),'VariableNames',{'Parameter','Original','Counterfactual'});
write_table(parameters,folder,'TableCF01_parameters','Production counterfactual and unchanged inputs');
baseline_rows=cell(0,9);policy_rows=cell(0,14);decomp_rows=cell(0,8);path_rows=cell(0,9);
for j=1:2
 model=models{j};b=model.baseline_cases{1};g=transition_abgp_rates(b.p,'paper');
 objects=resource_objects(b);lab=b.path.l_s(:)./(b.path.l_m(:)+b.path.l_s(:));
 baseline_rows(end+1,:)={labels{j},b.path.K(1)/objects.Y(1),lab(1),lab(15)-lab(1), ...
  objects.beta_bar(1),100*b.path.R_input(1)/b.path.R(1),100*b.path.R(end)/b.path.R(1),g.g_e,g.nonhom_growth_paper}; 
 for t=0:b.T
  path_rows(end+1,:)={labels{j},t,objects.beta_bar(t+1),objects.investment_share(t+1), ...
   lab(t+1),b.path.R(t+1)/b.path.R(1),objects.Y(t+1),b.path.R_input(t+1),b.path.h(t+1)}; 
 end
 for k=1:numel(model.policy_cases)
  r=model.policy_cases{k}.results{1};d=r.decomposition;
  partner=models{1}.policy_cases{find(cellfun(@(q)strcmp(q.id,r.case_name),models{1}.policy_cases))}.results{1};
  assert(isequaln(r.policy_specification,partner.policy_specification),'CF:Policy','Policy mismatch');
  cumulative=100*(cumsum(r.path.R_input(1:m.T))./cumsum(b.path.R_input(1:m.T))-1);
  H=m.cumulative_periods;v=zeros(1,numel(H));
  for h=1:numel(H)
   v(h)=cumulative(H(h));delta=sum(r.path.R_input(1:H(h)))-sum(b.path.R_input(1:H(h)));
   stock=b.path.R(H(h)+1)-r.path.R(H(h)+1);
   assert(abs(delta-stock)/max(1,abs(stock))<1e-5,'CF:StockFlow','Cumulative stock-flow check failed');
  end
  for h=1:numel(H)
   policy_rows(end+1,:)={r.case_name,labels{j},H(h),v(h), ...
    100*(r.path.R_input(1)/b.path.R_input(1)-1),r.effects.resource_intensity_effect_percent(1), ...
    d.total_log_points,d.component_log_points(1),d.component_log_points(2),d.component_log_points(3), ...
    d.dates(1),d.dates(2),sum(cumulative>1e-6),sum(cumulative< -1e-6)}; 
  end
  decomp_rows(end+1,:)={r.case_name,labels{j},d.dates(1),d.dates(2), ...
   d.component_log_points(1),d.component_log_points(2),d.component_log_points(3),d.total_log_points}; 
 end
end
baselines=cell2table(baseline_rows,'VariableNames',{'Scenario','CapitalOutputT0','ConditionalServiceLabourT0', ...
 'LabourShareChange0to14','BetaBarT0','InitialExtractionStockPercent','TerminalRemainingStockPercent','ABGPge','NonhomGrowth'});
write_table(baselines,folder,'TableCF02_baselines','Baseline changes without recalibration');
policies=cell2table(policy_rows,'VariableNames',{'Policy','Scenario','Periods','CumulativePercent', ...
 'InitialExtractionPercent','InitialIntensityPercent','TiltLogPoints','HotellingLogPoints','GrowthLogPoints', ...
 'StructuralLogPoints','EndpointStart','EndpointEnd','PositiveCumulativeWindows','NegativeCumulativeWindows'});
writetable(policies,fullfile(folder,'policy_comparison_full.csv'));
compact_rows=cell(0,9);
for k=1:numel(report.reference.policy_cases)
 id=report.reference.policy_cases{k}.id;
 for j=1:2
  mask=strcmp(policies.Policy,id)&strcmp(policies.Scenario,labels{j});z=policies(mask,:);
  get=@(H)z.CumulativePercent(z.Periods==H);
  if all(ismember([10,50,100,200],m.cumulative_periods))
   compact_rows(end+1,:)={id,labels{j},z.InitialExtractionPercent(1), ...
    z.InitialIntensityPercent(1),get(10),get(50),get(100),get(200),z.TiltLogPoints(1)}; 
  end
 end
end
if ~isempty(compact_rows)
 compact=cell2table(compact_rows,'VariableNames',{'Policy','Scenario','InitialExtractionPercent', ...
  'InitialIntensityPercent','Cumulative10Percent','Cumulative50Percent','Cumulative100Percent','Cumulative200Percent','TiltLogPoints'});
 write_table(compact,folder,'TableCF03_policy_effects','Policy effects relative to each model baseline');
end
decomposition=cell2table(decomp_rows,'VariableNames',{'Policy','Scenario','Date0','Date1','Hotelling','Growth','Structural','Total'});
write_table(decomposition,folder,'TableCF04_decomposition','Exact equilibrium decomposition in log points');
writetable(cell2table(path_rows,'VariableNames',{'Scenario','Date','BetaBar','InvestmentRevenueShare', ...
 'ConditionalServiceLabourShare','RemainingStockFraction','FinalOutput','Extraction','ResourceRent'}), ...
 fullfile(folder,'baseline_comparison_paths.csv'));
writetable(struct2table(report.channel_checks),fullfile(folder,'structural_channel_checks.csv'));
writetable(struct2table(report.cumulative_horizon_checks),fullfile(folder,'cumulative_horizon_checks.csv'));
horizon_rows=cell(0,7);
bc=report.flat.baseline_comparisons{1};
horizon_rows(end+1,:)={'baseline',bc.T,bc.reference_T,bc.screen_passed, ...
 100*bc.early_path.maximum_over_variables,NaN,NaN};
for k=1:numel(report.flat.policy_cases)
 r=report.flat.policy_cases{k};hc=r.comparisons{1};
 horizon_rows(end+1,:)={r.id,hc.T,hc.reference_T,hc.screen_passed, ...
  hc.early_maximum_real_effect_difference_pp,hc.extraction_t0_effect_difference_pp, ...
  max(abs(hc.decomposition_component_difference_log_points))}; 
end
writetable(cell2table(horizon_rows,'VariableNames',{'Case','T','CheckT','ScreenPassed', ...
 'EarlyDifference','InitialExtractionDifferencePP','DecompositionDifferenceLogPoints'}),fullfile(folder,'horizon_screen.csv'));
for j=1:2
 target=fullfile(out,strrep(lower(labels{j}),' ','_'));
 if ~isfolder(target),mkdir(target);end
 export_results(models{j},target);
end

fid=fopen(fullfile(folder,'counterfactual_parameters.m'),'w');assert(fid>=0,'CF:Export','Cannot save parameters');
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%% Complete counterfactual parameter choices, not a recalibration.\ncfg=solver_config();\n');
for k=1:9,fprintf(fid,'cfg.p.%s=%.17g;\n',keys{k},values(k,2));end
fprintf(fid,'cfg.p.zeta_s=-cfg.p.zeta_m;\ncfg.targets.K0=%.17g; cfg.targets.R0=%.17g;\n',values(10,2),values(11,2));
fprintf(fid,'%% Other preferences, A0 and growth remain those in solver_config.m.\n');clear cleanup
exports=struct('directory',folder,'parameters',parameters,'baselines',baselines,'policies',policies,'decomposition',decomposition);
save(fullfile(folder,'counterfactual_tables.mat'),'exports');

end

function o=resource_objects(r)
V=[r.path.p_m(:).*r.path.c_m(:),r.path.p_s(:).*r.path.c_s(:),r.path.I(:)];
o=struct('Y',sum(V,2),'investment_share',V(:,3)./sum(V,2), ...
 'beta_bar',V*[r.p.beta_m;r.p.beta_s;r.p.beta_x]./sum(V,2));
end
function r=get_policy(model,id)
k=find(cellfun(@(q)strcmp(q.id,id),model.policy_cases));assert(numel(k)==1,'CF:Plot','Policy missing');r=model.policy_cases{k}.results{1};
end
function write_table(t,folder,stem,caption)
writetable(t,fullfile(folder,[stem '.csv']));
fid=fopen(fullfile(folder,[stem '.tex']),'w');assert(fid>=0,'CF:Export','Cannot save table');cleanup=onCleanup(@()fclose(fid)); 
fprintf(fid,'\\begin{table}[htbp]\n\\centering\n\\caption{%s}\n\\label{tab:%s}\n',caption,strrep(stem,'_','-'));
fprintf(fid,'\\resizebox{\\textwidth}{!}{%%\n\\begin{tabular}{%s}\n\\toprule\n',repmat('l',1,width(t)));
headers=regexprep(t.Properties.VariableNames,'([a-z])([A-Z])','$1 $2');headers=strrep(headers,'_','\_');
fprintf(fid,'%s \\\\\n\\midrule\n',strjoin(headers,' & '));
for i=1:height(t)
 row=cell(1,width(t));
 for j=1:width(t)
  v=t{i,j};if iscell(v),v=v{1};end
  if isnumeric(v)||islogical(v)
   if ~isfinite(v),row{j}='--';elseif abs(v)<1e-10,row{j}='0';
   elseif abs(v)<1e-3||abs(v)>=1e4,row{j}=sprintf('%.4e',v);else,row{j}=sprintf('%.6f',v);end
  else,row{j}=strrep(strrep(char(v),'_','\_'),'&','\&');end
 end
 fprintf(fid,'%s \\\\\n',strjoin(row,' & '));
end
fprintf(fid,'\\bottomrule\n\\end{tabular}}\n\\end{table}\n');
end
