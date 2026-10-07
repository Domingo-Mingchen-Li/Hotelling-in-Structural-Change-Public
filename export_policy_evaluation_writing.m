function summary=export_policy_evaluation_writing(outputs_root,options)
if nargin<1||isempty(outputs_root),outputs_root=fullfile(fileparts(mfilename('fullpath')),'outputs');end
if nargin<2||isempty(options),options=struct();end
options=defaults(options,struct('destination','','T',200,'eta',2,'loss_target',.01, ...
 'last_local_date',80,'make_plots',true,'make_tables',true,'vector_pdf',true,'resolution',300,'line_width',1.6, ...
 'policy_ids',{{}},'climate_report_file','', ...
 'welfare_report_file','','ramsey_report_file','','shared_report_file',''));
o=options;out=new_folder(outputs_root,'writing_policy_evaluation',o);
[c,cf]=accepted_report(fullfile(outputs_root,'climate_parameterization'),'climate_report.mat',o.climate_report_file);
[w,wf]=accepted_report(fullfile(outputs_root,'industrial_welfare'),'industrial_welfare_report.mat',o.welfare_report_file);
[r,rf]=accepted_report(fullfile(outputs_root,'ramsey'),'ramsey_report.mat',o.ramsey_report_file);
[s,sf]=accepted_report(fullfile(outputs_root,'ramsey_shared_state_counterfactual'),'ramsey_shared_state_report.mat',o.shared_report_file);
assert(isequaln(w.common_inherited,r.common_inherited)&&isequaln(s.counterfactual.common_inherited,r.common_inherited), ...
 'Writing:States','Welfare/Ramsey/shared-state reports have different inherited states');
assert(isequaln(w.common_normalization,r.common_normalization)&&isequaln(c.main.anchors,r.common_normalization), ...
 'Writing:Normalization','Reports use different environmental normalization');
assert(r.scenario.eta==o.eta&&abs(r.scenario.loss_target-o.loss_target)<1e-12, ...
 'Writing:Scenario','Selected Ramsey scenario missing; select matching existing reports');
ix=find([c.main.scenarios.eta]==o.eta&abs([c.main.scenarios.loss_target]-o.loss_target)<1e-12,1);
assert(~isempty(ix)&&abs(c.main.scenarios(ix).chi/r.scenario.chi-1)<1e-10,'Writing:Chi','Climate/Ramsey chi mismatch');
assert(isequaln(s.reference.common_inherited,r.common_inherited)&&abs(s.reference.scenario.chi/r.scenario.chi-1)<1e-12, ...
 'Writing:Reference','Shared-state reference differs');
sources={cf;wf;rf;sf};types={'Climate parameters';'Industrial welfare';'Original Ramsey';'Shared-state Ramsey counterfactual'};
m=c.configuration;
kernel=table((1:numel(m.a))',m.a(:),m.timescales_years(:),m.rho(:), ...
 'VariableNames',{'Component','Weight','TimescaleYears','Persistence'});
if o.make_tables,paper_table(kernel,out,'TableA4','Reduced climate kernel');end
parameters=struct2table(c.main.scenarios);
if o.make_tables,paper_table(parameters(:,{'eta','loss_target','utility_loss','chi'}),out, ...
 'TableA5','Normative environmental damage scenarios');end
anchors=table(m.L,c.main.anchors.e_reference,c.main.anchors.M_reference,c.main.anchors.reference_date, ...
 c.main.benchmark.utility_loss,'VariableNames',{'AnnouncementDate','ResourceNormalization','BurdenNormalization','BurdenReferenceDate','FixedUtilityLossAnchor'});

wi=find(cellfun(@(x)x.calendar_terminal==o.T,w.runs),1);assert(~isempty(wi),'Writing:Horizon','Welfare horizon missing');
allrows=struct2table(w.runs{wi}.rows);

selected=allrows(allrows.eta==o.eta&abs(allrows.loss_target-o.loss_target)<1e-12,:);
assert(~isempty(selected),'Writing:Scenario','No selected welfare scenario');
if ~isempty(o.policy_ids)
 assert(all(ismember(string(o.policy_ids),string(selected.policy))),'Writing:Policy','Unknown policy ID');
 selected=selected(ismember(string(selected.policy),string(o.policy_ids)),:);
end
assert(max(abs(selected.chi/r.scenario.chi-1))<1e-10,'Writing:Chi','Welfare chi differs');
main=selected(:,{'policy','consumption_utility_change','environment_welfare_change','net_welfare_change', ...
 'consumption_in_anchor_units','environment_in_anchor_units','net_in_anchor_units','finite_cumulative_extraction_percent'});
if o.make_tables,paper_table(main(:,{'policy','consumption_in_anchor_units','environment_in_anchor_units','net_in_anchor_units'}),out,'TableA6','Industrial policy welfare effects');end

if o.make_plots
 f=paper_figure(1,1);set(f,'Position',[2 2 23 12]);cleanup=onCleanup(@()close(f));
 bars=bar([selected.consumption_in_anchor_units selected.environment_in_anchor_units selected.net_in_anchor_units]);welfare_colors(bars);
 yline(0,'k:','HandleVisibility','off');grid on;xticks(1:height(selected));xticklabels(strrep(cellstr(string(selected.policy)),'_',' '));xtickangle(20);
 ylabel('$\Delta W_j/\mathcal{A}$');paper_legend(bars,{'Consumption','Environment','Net'});
 title(sprintf('Industrial-policy welfare: calendar T=%d',o.T));
 save_plot(f,out,'Figure14',o);clear cleanup
end
rt=readtable(fullfile(fileparts(rf),'ramsey_summary.csv'));

p=readtable(fullfile(fileparts(rf),sprintf('ramsey_path_T%d.csv',o.T)));

H=height(p)-1;assert(o.last_local_date>=1&&fix(o.last_local_date)==o.last_local_date,'Writing:Dates','Invalid local plot window');
cum=100*(cumsum(p.policy_extraction(1:H))./cumsum(p.baseline_extraction(1:H))-1);
periods=unique([10 50 100 H]);periods=periods(periods<=H);
flows=table(periods(:),p.calendar_date(periods(:))+1,cum(periods(:)), ...
 'VariableNames',{'PostAnnouncementFlowPeriods','CalendarEndpoint','CumulativeExtractionPercent'});

if o.make_plots
 f=paper_figure(2,2);cleanup=onCleanup(@()close(f));tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
 D=min(o.last_local_date+1,height(p));
 nexttile;tax_lines=plot(p.calendar_date(1:D),100*[p.resource_tax(1:D) p.conditional_pigou_tax(1:D)],'LineWidth',o.line_width);
 ylabel('$100\tau_{e,t}$');xlabel('$t$');paper_legend(tax_lines,{'Numerical optimum','Pigouvian diagnostic'});grid on;title('Committed resource tax');
 nexttile;plot(p.calendar_date,100*(p.policy_extraction./p.baseline_extraction-1),'LineWidth',o.line_width);yline(0,'k:','HandleVisibility','off');ylabel('$100(e_t/e_t^b-1)$');xlabel('$t$');grid on;title('Extraction relative to baseline');
 nexttile;plot(p.calendar_date,100*(p.policy_burden./p.baseline_burden-1),'LineWidth',o.line_width);yline(0,'k:','HandleVisibility','off');ylabel('$100(M_t/M_t^b-1)$');xlabel('$t$');grid on;title('Environmental burden relative to baseline');
 nexttile;plot(1:H,cum,'LineWidth',o.line_width);yline(0,'k:','HandleVisibility','off');xlabel('$H$');ylabel('$100(\mathcal{E}_H/\mathcal{E}_H^b-1)$');grid on;title('Cumulative extraction');
 save_plot(f,out,'Figure15',o);clear cleanup
end
[shared,sharedpath]=comparison_exports(s,sf,out,o,'Shared','Common inherited states');
manifest=table(string(types(:)),string(sources(:)),'VariableNames',{'Content','SourceReport'});

summary=struct('output_directory',out,'options',o,'sources',manifest,'industrial_welfare',main, ...
 'ramsey_welfare',rt,'shared_state_comparison',shared,'shared_state_paths',sharedpath);

end
function [t,path]=comparison_exports(r,file,out,o,tag,scope)
assert(r.passed,'Writing:Acceptance','Counterfactual report failed');
if strcmp(tag,'Shared')
 assert(r.shared_state_audit.passed&&isequaln(r.counterfactual.common_inherited,r.reference.common_inherited), ...
  'Writing:States','Shared-state equality check failed');
end
folder=fileparts(file);t=readtable(fullfile(folder,'optimal_welfare_comparison.csv'));
t.InitialTaxPercent=100*t.initial_resource_tax;
if o.make_tables,paper_table(t(:,{'economy','calendar_terminal','InitialTaxPercent','consumption_utility_change', ...
 'environment_welfare_change','net_welfare_change'}),out,'TableA7',['Optimal policy comparison: ' scope]);end
path=readtable(fullfile(folder,sprintf('optimal_tax_comparison_T%d.csv',o.T)));

inherited=readtable(fullfile(folder,'inherited_states_comparison.csv'));

if ~o.make_plots,return;end
f=paper_figure(2,2);cleanup=onCleanup(@()close(f)); 
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');D=min(o.last_local_date+1,height(path));dates=path.calendar_date(1:D);
nexttile;tax_lines=plot(dates,100*[path.original_optimal_resource_tax(1:D),path.flat_optimal_resource_tax(1:D)],'LineWidth',o.line_width);
paper_legend(tax_lines,{'Original','Flat input shares'});ylabel('$100\tau_{e,t}$');xlabel('$t$');title(scope);grid on;
nexttile;plot(dates,[path.original_extraction_deviation_percent(1:D),path.flat_extraction_deviation_percent(1:D)],'LineWidth',o.line_width);
yline(0,'k:','HandleVisibility','off');ylabel('$100(e_t/e_t^b-1)$');xlabel('$t$');title('Extraction response');grid on;
nexttile;ingredient_lines=plot(dates,[path.flat_marginal_environment_cost(1:D)./path.original_marginal_environment_cost(1:D), ...
 path.flat_mu_times_h(1:D)./path.original_mu_times_h(1:D)],'LineWidth',o.line_width);yline(1,'k:','HandleVisibility','off');
paper_legend(ingredient_lines,{'$\mathrm{MEC}_t$','$\mu_t h_t$'},'latex');ylabel('$X_t^{\mathrm{flat}}/X_t^{\mathrm{orig}}$');xlabel('$t$');title('Conditional Pigouvian ingredients');grid on;
z=t(t.calendar_terminal==o.T,:);assert(height(z)==2,'Writing:Comparison','Expected two economies');
nexttile;bars=bar([z.consumption_utility_change z.environment_welfare_change z.net_welfare_change]);welfare_colors(bars);
xticklabels({'Original','Flat input shares'});paper_legend(bars,{'Consumption','Environment','Net'});ylabel('$\Delta W_j$');yline(0,'k:','HandleVisibility','off');grid on;
save_plot(f,out,'Figure16',o);
end
function [r,file]=accepted_report(root,name,explicit)
if ~isempty(explicit)
 file=explicit;assert(isfile(file),'Writing:Source','Missing explicit report');z=load(file,'report');r=z.report;
 assert(r.passed,'Writing:Acceptance','Explicit report has not passed');return;
end
files=dir(fullfile(root,'**',name));accepted=[];
for i=1:numel(files)
 z=load(fullfile(files(i).folder,files(i).name),'report');
 if isfield(z,'report')&&isfield(z.report,'passed')&&z.report.passed,accepted(end+1)=i;end 
end
assert(~isempty(accepted),'Writing:Source','No accepted %s under %s',name,root);
[~,j]=max([files(accepted).datenum]);k=accepted(j);file=fullfile(files(k).folder,files(k).name);z=load(file,'report');r=z.report;

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
function lg=paper_legend(handles,labels,interpreter)
if nargin<3,interpreter='none';end
lg=legend(handles,labels,'AutoUpdate','off','Location','southoutside', ...
 'Orientation','horizontal','NumColumns',numel(labels), ...
 'Interpreter',interpreter,'Box','off','FontName','Times New Roman','FontSize',9);
lg.ItemTokenSize=[14 10];
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

function welfare_colors(bars)
colors=[.08 .28 .48;.20 .48 .31;.47 .32 .60];
for k=1:numel(bars),bars(k).FaceColor=colors(k,:);end
end
