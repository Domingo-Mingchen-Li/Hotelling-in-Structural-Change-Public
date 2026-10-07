function summary = export_measurement_writing(outputs_root, options)
if nargin<1||isempty(outputs_root),outputs_root=fullfile(fileparts(mfilename('fullpath')),'outputs');end
if nargin<2||isempty(options),options=struct();end
options=defaults(options,struct('destination','','case_ids',{{'annual_average','domestic_inputs','mapping_c26'}}, ...
 'case_labels',{{'Annual average','Domestic inputs','C26 mapping'}}, ...
 'last_date',40,'make_plots',true,'resolution',300,'line_width',1.6,'vector_pdf',true));
assert(numel(options.case_ids)==numel(options.case_labels),'Writing:Cases','IDs and labels must match');
assert(options.last_date>=1&&options.last_date<=200&&fix(options.last_date)==options.last_date,'Writing:Dates','last_date must be an integer from 1 to 200');
destination=new_folder(outputs_root,'writing_measurement',options);
ids=options.case_ids;labels=options.case_labels;case_status=cell(0,3);
allfiles=dir(fullfile(outputs_root,'**','measurement_report.mat'));
records={};names={};source={};alternative={};reference=[];
for k=1:numel(ids)
 candidates={};dates=[];
 for j=1:numel(allfiles)
  file=fullfile(allfiles(j).folder,allfiles(j).name);raw=load(file,'report');
  if ~isfield(raw,'report')||~isfield(raw.report,'case_id')||~strcmp(raw.report.case_id,ids{k}),continue;end
  q=raw.report;
  if ~q.passed||~q.calibration_report.passed,continue;end
  candidates{end+1}=file;dates(end+1)=allfiles(j).datenum; 
 end
 if isempty(candidates)
  case_status(end+1,:)={ids{k},false,'No accepted saved report'};continue;
 end
 [~,latest]=max(dates);file=candidates{latest};raw=load(file,'report');q=raw.report;
 assert(q.calibration_report.fit_passed&&q.calibration_report.identification_passed&& ...
  q.calibration_report.finite_horizon_passed&&q.calibration_report.horizon_passed, ...
  'E:Audit','Incomplete calibration acceptance: %s',file);
 assert(q.configuration.T==200,'E:Horizon','This exporter expects T=200');
 alt=load_solve(fullfile(fileparts(file),'report.mat'));
 refdir=resolve_reference(outputs_root,q.reference_output_dir);
 ref=load_solve(fullfile(refdir,'report.mat'));
 assert(ref.passed&&alt.passed,'E:Audit','Policy results did not pass');
 fitted=[alt.configuration.p.zeta_m,alt.configuration.targets.K0,alt.configuration.targets.R0];
 assert(max(abs(fitted-q.calibration_report.parameters(:).')./max(1,abs(fitted)))<1e-8, ...
  'E:Fit','Policy report does not use the reported fitted parameters');
 for sector={'m','s','x','e'}
  for kind={'alpha','beta'}
   key=[kind{1} '_' sector{1}];
   assert(alt.configuration.p.(key)==q.calibration_report.configuration.base.p.(key), ...
    'E:Fit','Production inputs differ between fit and policies');
  end
 end
 if isempty(reference),reference=ref;else
  assert(isequaln(reference.configuration.p,ref.configuration.p)&& ...
   isequaln(reference.configuration.targets,ref.configuration.targets)&& ...
   isequaln(reference.configuration.A0,ref.configuration.A0), ...
   'E:Reference','Cases use different paper references; export them separately');
 end
 compare=struct2table(q.comparison.rows);
 for j=1:height(compare)
  metric=char(compare.metric(j));pid=char(compare.policy(j));
  if ismember(metric,{'total_log_points','structural_log_points','cumulative_10_periods_percent'})
   refvalue=saved_metric(ref,pid,metric);altvalue=saved_metric(alt,pid,metric);
   assert(abs(altvalue-compare.alternative(j))<1e-7*max(1,abs(altvalue)), ...
    'E:Fit','Saved comparison disagrees with alternative for %s',pid);
   assert(abs(refvalue-compare.reference(j))<1e-7*max(1,abs(refvalue)), ...
    'E:Reference','Saved comparison disagrees with reference for %s',pid);
  end
 end
 records{end+1}=q;alternative{end+1}=alt;names{end+1}=labels{k};source{end+1}=file; 
 case_status(end+1,:)={ids{k},true,file};
end
assert(~isempty(records),'E:Missing','No accepted measurement results found');
raw=load(fullfile(outputs_root,'calibration_output','calibration_report.mat'),'report');original=raw.report;
assert(original.passed,'E:Reference','Original calibration did not pass');
expected=[reference.configuration.p.zeta_m,reference.configuration.targets.K0,reference.configuration.targets.R0];
assert(max(abs(original.parameters(:).'-expected)./max(1,abs(expected)))<1e-8, ...
 'E:Reference','Original calibration parameters disagree with policy reference');
fits=[{original},cellfun(@(q)q.calibration_report,records,'UniformOutput',false)];
models=[{reference},alternative];case_names=[{'Paper baseline'},names];
parameter_keys={'alpha_m','beta_m','alpha_s','beta_s','alpha_x','beta_x','alpha_e','beta_e','zeta_m','K0','R0'};
values=zeros(numel(parameter_keys),numel(models));
for j=1:numel(models)
 cfg=models{j}.configuration;
 for k=1:9,values(k,j)=cfg.p.(parameter_keys{k});end
 values(10:11,j)=[cfg.targets.K0;cfg.targets.R0];
end
parameters=array2table(values,'VariableNames',matlab.lang.makeValidName(case_names));
parameters=addvars(parameters,string(parameter_keys(:)),'Before',1,'NewVariableNames','Parameter');
write_outputs(parameters,destination,'TableA2','Parameters across measurement specifications');
moment_rows=cell(0,5);audit_rows=cell(0,15);
for j=1:numel(fits)
 c=fits{j};base=models{j}.baseline_cases{1};e=base.path.R_input(:);R=base.path.R(:);
 for k=1:3
  moment_rows(end+1,:)={case_names{j},c.configuration.moment_names{k}, ...
   c.configuration.targets(k),c.moments(k),c.moment_errors(k)}; 
 end
 audit_rows(end+1,:)={case_names{j},c.passed,c.fit_passed,c.identification_passed, ...
  c.finite_horizon_passed,c.horizon_passed,c.identification.rank, ...
  c.identification.condition_number,min(c.identification.singular_values), ...
  c.multistart.maximum_spread,max(abs(c.horizon.difference)), ...
  100*e(1)/R(1),100*R(end)/R(1),c.configuration.lower(2),c.configuration.lower(3)}; 
end
moments=cell2table(moment_rows,'VariableNames',{'Specification','Moment','Target','Model','Error'});

audits=cell2table(audit_rows,'VariableNames',{'Specification','OverallPassed','FitPassed', ...
 'IdentificationPassed','EquilibriumPassed','HorizonPassed','JacobianRank','JacobianCondition', ...
 'SmallestSingularValue','MultistartSpread','MaxHorizonMomentDifference', ...
 'InitialExtractionStockPercent','TerminalRemainingStockPercent','LowerK0','LowerR0'});

pids=cellfun(@(q)q.id,reference.policy_cases,'UniformOutput',false);
policy_rows=cell(0,7);decomp_rows=cell(0,7);
for j=1:numel(models)
 for k=1:numel(pids)
  pid=pids{k};r=policy_result(models{j},pid);b=models{j}.baseline_cases{1};
  reference_policy=policy_result(reference,pid);
  assert(isequaln(r.policy_specification,reference_policy.policy_specification), ...
   'E:Policy','Policy schedules differ across specifications');
  cumulative=100*(cumsum(r.path.R_input(1:200))./cumsum(b.path.R_input(1:200))-1);
  policy_rows(end+1,:)={policy_label(pid),case_names{j},cumulative(10), ...
   cumulative(50),cumulative(100),cumulative(200),r.effects.resource_intensity_effect_percent(1)}; 
  d=r.decomposition;
  assert(abs(sum(d.component_log_points)-d.total_log_points)<1e-7,'E:Decomposition','Identity failed');
  if ~strcmp(pid,'consumption_permanent')
   decomp_rows(end+1,:)={policy_label(pid),case_names{j},d.component_log_points(1), ...
    d.component_log_points(2),d.component_log_points(3),d.total_log_points, ...
    100*d.signed_component_shares(3)}; 
  end
 end
end
policies=cell2table(policy_rows,'VariableNames',{'Policy','Specification','Cumulative10Percent', ...
 'Cumulative50Percent','Cumulative100Percent','Cumulative200Percent','InitialIntensityPercent'});
values=zeros(numel(pids),numel(case_names));
for j=1:numel(case_names)
 for k=1:numel(pids)
  rows=strcmp(string(policies.Policy),policy_label(pids{k}))&strcmp(string(policies.Specification),case_names{j});
  values(k,j)=policies.Cumulative10Percent(rows);
 end
end
t3=array2table(values,'VariableNames',matlab.lang.makeValidName(case_names));
labels=string(cellfun(@policy_label,pids,'UniformOutput',false));
t3=addvars(t3,labels(:),'Before',1,'NewVariableNames','Policy');
write_outputs(t3,destination,'TableA3','Cumulative extraction over the first ten periods');
decomposition=cell2table(decomp_rows,'VariableNames',{'Policy','Specification','HotellingLogPoints', ...
 'GrowthLogPoints','StructuralLogPoints','TotalLogPoints','StructuralSignedSharePercent'});

plot_ids=pids(~strcmp(pids,'consumption_permanent'));
assert(numel(plot_ids)<=6,'E:Plot','Exporter supports at most six nonneutral policies');
if options.make_plots
for metric=1:2
 f=paper_figure(2,3);
 cleanup=onCleanup(@()close(f));
 layout=tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');
 palette=[.08 .28 .48;.74 .27 .16;.20 .48 .31;.47 .32 .60;.45 .45 .45];
 styles={'-','--',':','-.','--'};
 assert(numel(models)<=size(palette,1),'Writing:Palette','Extend palette for additional cases');
 for k=1:numel(plot_ids)
  ax=nexttile(layout);hold(ax,'on');
  for j=1:numel(models)
   r=policy_result(models{j},plot_ids{k});b=models{j}.baseline_cases{1};
   if metric==1,y=100*(r.path.R_input(:)./b.path.R_input(:)-1);
   else,y=r.effects.resource_intensity_effect_percent(:);end
   plot(ax,0:options.last_date,y(1:options.last_date+1),'Color',palette(j,:),'LineStyle',styles{j},'LineWidth',options.line_width, ...
    'DisplayName',case_names{j});
  end
  yline(ax,0,':','Color',[.65,.65,.65],'HandleVisibility','off');
  if ~r.policy_specification.permanent
   xline(ax,9.5,':','Color',[.65,.65,.65],'HandleVisibility','off');
  end
  title(ax,policy_label(plot_ids{k}));xlabel(ax,'$t$');if metric==1,ylabel(ax,'$100(e_t/e_t^b-1)$');else,ylabel(ax,'$100[(e_t/Y_t)/(e_t^b/Y_t^b)-1]$');end
  grid(ax,'on');xlim(ax,[0,options.last_date]);set(ax,'FontSize',10);
  if k==1,legend(ax,'Location','best','FontSize',8);end
 end
 if metric==1,stem='FigureA3';title(layout,'Resource extraction relative to each specification baseline');
 else,stem='FigureA4';title(layout,'Resource intensity relative to each specification baseline');end
 save_plot(f,destination,stem,options);
 clear cleanup
end
end

summary=struct('output_directory',destination,'options',options,'included_cases',{names},'source_reports',{source}, ...
 'parameters',parameters,'moments',moments,'validation',audits,'policies',policies, ...
 'decomposition',decomposition,'cumulative_convention','H flows at t=0,...,H-1; terminal flow excluded', ...
 'decomposition_convention','Log points; equilibrium accounting, not isolated causal effects');

end

function r=load_solve(file)
assert(isfile(file),'E:File','Saved solver report missing: %s',file);
s=load(file,'report');r=s.report;
assert(strcmp(r.model_mode,'paper')&&isequal(r.horizons,200),'E:Model','Expected paper T=200 results');
end
function folder=resolve_reference(root,saved)
parts=regexp(char(saved),'[\\/]','split');parts=parts(~cellfun('isempty',parts));name=parts{end};
d=dir(fullfile(root,'**','report.mat'));matches={};
for k=1:numel(d)
 [~,basename]=fileparts(d(k).folder);
 if strcmp(basename,name),matches{end+1}=d(k).folder;end 
end
assert(numel(matches)==1,'E:Reference','Expected one saved reference run %s; found %d',name,numel(matches));
folder=matches{1};
end
function r=policy_result(model,id)
k=find(cellfun(@(q)strcmp(q.id,id),model.policy_cases));
assert(numel(k)==1,'E:Policy','Missing or duplicate policy %s',id);
r=model.policy_cases{k}.results{1};assert(r.passed,'E:Policy','Policy not validated: %s',id);
end
function value=saved_metric(model,id,metric)
r=policy_result(model,id);b=model.baseline_cases{1};
switch metric
 case 'total_log_points',value=r.decomposition.total_log_points;
 case 'structural_log_points',value=r.decomposition.component_log_points(3);
 case 'cumulative_10_periods_percent',value=100*(sum(r.path.R_input(1:10))/sum(b.path.R_input(1:10))-1);
 otherwise,error('E:Metric','Unknown metric');
end
end
function label=policy_label(id)
switch id
 case 'service_temporary',label='Temporary services subsidy';
 case 'service_permanent',label='Permanent services subsidy';
 case 'consumption_temporary',label='Temporary consumption tax';
 case 'consumption_permanent',label='Permanent consumption tax';
 case 'carbon_temporary',label='Temporary carbon tax';
 case 'investment_purchase_20',label='Investment subsidy (20 percent)';
 case 'investment_purchase_matched',label='Matched investment subsidy';
 otherwise,label=strrep(id,'_',' ');
end
end
function write_outputs(t,folder,stem,caption)
writetable(t,fullfile(folder,[stem '.csv']));
fid=fopen(fullfile(folder,[stem '.tex']),'w');assert(fid>=0,'E:Output','Cannot save table');
cleanup=onCleanup(@()fclose(fid)); 
fprintf(fid,'\\begin{table}[htbp]\n\\centering\n\\caption{%s}\n\\label{tab:%s}\n',caption,strrep(stem,'_','-'));
fprintf(fid,'\\resizebox{\\textwidth}{!}{%%\n\\begin{tabular}{%s}\n\\toprule\n',repmat('l',1,width(t)));
head=regexprep(t.Properties.VariableNames,'([a-z])([A-Z])','$1 $2');
fprintf(fid,'%s \\\\\n\\midrule\n',strjoin(cellfun(@tex_escape,head,'UniformOutput',false),' & '));
for i=1:height(t)
 row=cell(1,width(t));
 for j=1:width(t)
  v=t{i,j};if iscell(v),v=v{1};end
  if isnumeric(v)||islogical(v)
   if ~isfinite(v),row{j}='--';
   elseif v==0||abs(v)<1e-10,row{j}='0';
   elseif abs(v)<1e-3||abs(v)>=1e4,row{j}=sprintf('%.4e',v);
   else,row{j}=sprintf('%.6f',v);end
  else,row{j}=tex_escape(char(v));end
 end
 fprintf(fid,'%s \\\\\n',strjoin(row,' & '));
end
fprintf(fid,'\\bottomrule\n\\end{tabular}}\n\\end{table}\n');
end
function s=tex_escape(s)
s=strrep(s,'_','\_');s=strrep(s,'%','\%');s=strrep(s,'&','\&');
end

function o=defaults(o,d)
for f=fieldnames(d).',if ~isfield(o,f{1})||isempty(o.(f{1})),o.(f{1})=d.(f{1});end,end
end
function folder=new_folder(root,tag,o)
if ~isfolder(root),error('Writing:Root','Outputs directory missing: %s',root);end
folder=o.destination;
if isempty(folder)
 d=dir(fullfile(root,[tag '_*']));d=d([d.isdir]);
 if isempty(d),folder=fullfile(root,[tag '_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);
 else,[~,k]=max([d.datenum]);folder=fullfile(d(k).folder,d(k).name);end
end
if ~isfolder(folder),mkdir(folder);end

end
function save_plot(f,out,stem,o)
paper_style(f,o);drawnow;
exportgraphics(f,fullfile(out,[stem '.png']),'Resolution',o.resolution,'BackgroundColor','white');
if o.vector_pdf,exportgraphics(f,fullfile(out,[stem '.pdf']),'ContentType','vector','BackgroundColor','white');end

end
function f=paper_figure(rows,columns)
f=figure('Visible','off','Color','w','Units','centimeters', ...
 'Position',[2 2 7*columns+2 5.5*rows+2], ...
 'DefaultAxesColorOrder',[.08 .28 .48;.74 .27 .16;.20 .48 .31;.47 .32 .60]);
end
function paper_style(f,o)
text_objects=findall(f,'Type','text');
for h=reshape(text_objects,1,[]),set(h,'FontName','Times New Roman','FontSize',10);end
palette=[.08 .28 .48;.74 .27 .16;.20 .48 .31;.47 .32 .60];
axes=findall(f,'Type','axes');
for ax=reshape(axes,1,[])
 set(ax,'FontName','Times New Roman','FontSize',10,'TickDir','out','Box','off');
 ax.YAxis.Exponent=0;grid(ax,'on');
 set(ax.Title,'FontWeight','normal','Interpreter','none');
 set(ax.XLabel,'Interpreter','latex');set(ax.YLabel,'Interpreter','latex');
 lines=findall(ax,'Type','line');
 for h=reshape(lines,1,[])
  if strcmp(h.HandleVisibility,'off'),continue;end
  h.LineWidth=o.line_width;
 end
end
legends=findall(f,'Type','legend');
for h=reshape(legends,1,[]),set(h,'Box','off','FontName','Times New Roman','FontSize',9);end
end
function paper_table(t,out,stem,caption)
writetable(t,fullfile(out,[stem '.csv']));
fid=fopen(fullfile(out,[stem '.tex']),'w');assert(fid>=0,'Writing:Output','Cannot write table');
c=onCleanup(@()fclose(fid)); 
fprintf(fid,'\\begin{table}[htbp]\n\\centering\n\\caption{%s}\n\\label{tab:%s}\n',escape_tex(caption),strrep(stem,'_','-'));
fprintf(fid,'\\resizebox{\\textwidth}{!}{%%\n\\begin{tabular}{%s}\n\\toprule\n',repmat('l',1,width(t)));
head=regexprep(t.Properties.VariableNames,'([a-z])([A-Z])','$1 $2');
fprintf(fid,'%s \\\\\n\\midrule\n',strjoin(cellfun(@escape_tex,head,'UniformOutput',false),' & '));
for i=1:height(t)
 row=cell(1,width(t));
 for j=1:width(t)
  v=t{i,j};if iscell(v),v=v{1};end
  if isnumeric(v)||islogical(v)
   if isnan(v),row{j}='--';elseif isinf(v),row{j}='$\infty$';
   elseif v==0,row{j}='0';elseif abs(v)<1e-4||abs(v)>=1e5,row{j}=sprintf('%.4e',v);
   else,row{j}=sprintf('%.6g',v);end
  else,row{j}=escape_tex(char(string(v)));end
 end
 fprintf(fid,'%s \\\\\n',strjoin(row,' & '));
end
fprintf(fid,'\\bottomrule\n\\end{tabular}}\n\\end{table}\n');
end
function s=escape_tex(s)
s=strrep(s,'_','\_');s=strrep(s,'%','\%');s=strrep(s,'&','\&');s=strrep(s,'#','\#');
end
