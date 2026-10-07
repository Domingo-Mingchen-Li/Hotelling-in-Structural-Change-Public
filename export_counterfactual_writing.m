function summary = export_counterfactual_writing(outputs_root,options)
if nargin<1||isempty(outputs_root),outputs_root=fullfile(fileparts(mfilename('fullpath')),'outputs');end
if nargin<2||isempty(options),options=struct();end
options=defaults(options,struct('destination','','report_file','','last_date',40, ...
 'policy_ids',{{'service_temporary','consumption_temporary'}},'make_plots',true,'resolution',300,'line_width',1.6,'vector_pdf',true));
source=fullfile(outputs_root,'resource_intensity_counterfactual');
if ~isempty(options.report_file),source=options.report_file;end
report_file=find_report(source);z=load(report_file,'report');report=z.report;
assert(report.passed&&report.channel_closed&&report.horizon_screen_passed&&report.cumulative_horizon_passed, ...
 'Writing:Acceptance','Counterfactual acceptance not complete');
T=report.configuration.T;last_date=options.last_date;
assert(T>=200&&last_date>=1&&last_date<=T&&fix(last_date)==last_date,'Writing:Dates','Invalid plotting horizon');
destination=new_folder(outputs_root,'writing_counterfactual',options);
models={report.reference,report.flat};names={'Original model','Equal resource shares'};
colors=[.08 .28 .48;.74 .27 .16];styles={'-','--'};
ids=cellfun(@(q)q.id,models{1}.policy_cases,'UniformOutput',false);
if ~isempty(options.policy_ids),assert(all(ismember(options.policy_ids,ids)),'Writing:Policy','Unknown policy');ids=options.policy_ids;end
data=cell(numel(ids),2);metrics=cell(0,16);path_rows=cell(0,12);
for k=1:numel(ids)
 for j=1:2
  r=get_case(models{j},ids{k},T);b=get_baseline(models{j},T);
  if j==2
   reference_case=get_case(models{1},ids{k},T);
   assert(isequaln(r.policy_specification,reference_case.policy_specification), ...
    'CFPlot:Policy','Different policy specifications');
  end
  u=policy_data(r,b,T);data{k,j}=u;d=r.decomposition;
  assert(max(abs(u.intensity-r.effects.resource_intensity_effect_percent(:)))<1e-7, ...
   'CFPlot:Intensity','Recomputed resource intensity differs from saved report');
  assert(abs(sum(d.component_log_points)-d.total_log_points)<1e-7,'CFPlot:Identity','Decomposition identity failed');
  if j==2,assert(max(abs(u.beta_bar/report.common_beta-1))<1e-12,'CFPlot:Shares','Flat shares are not constant');end
  metrics(end+1,:)={ids{k},names{j},u.extraction(1),u.intensity(1), ...
   u.cumulative(10),u.cumulative(50),u.cumulative(100),u.cumulative(T), ...
   d.total_log_points,d.component_log_points(1),d.component_log_points(2),d.component_log_points(3), ...
   d.dates(1),d.dates(2),u.initial_stock_flow_error,max(abs(u.gross_return))}; 
  for t=0:T
   ret=NaN;if t<T,ret=u.gross_return(t+1);end
   cumulative=NaN;if t>=1,cumulative=u.cumulative(t);end
   path_rows(end+1,:)={ids{k},names{j},t,u.extraction(t+1),u.intensity(t+1), ...
    u.K(t+1),u.I(t+1),u.investment_share_pp(t+1),u.service_labour_pp(t+1), ...
    ret,cumulative,u.beta_bar(t+1)}; 
  end
 end
 if options.make_plots
 f=paper_figure(3,3);
 cleanup=onCleanup(@()close(f));layout=tiledlayout(f,3,3,'TileSpacing','compact','Padding','compact');
 keys={'extraction','intensity','cumulative','K','I','investment_share_pp','service_labour_pp','gross_return'};
 titles={'Resource extraction','Resource intensity','Cumulative extraction', ...
  'Physical capital','Investment output','Investment share of final output', ...
  'Services share of m+s labour','Gross return to investment'};
 for n=1:8
  ax=nexttile(layout);hold(ax,'on');
  for j=1:2
   y=data{k,j}.(keys{n});
   if n==3,x=1:last_date;index=x;
   elseif n==8,x=0:min(last_date,T-1);index=x+1;
   else,x=0:last_date;index=x+1;end
   plot(ax,x,y(index),'Color',colors(j,:),'LineStyle',styles{j},'LineWidth',options.line_width,'DisplayName',names{j});
  end
  units=math_axis(keys{n});
  ylabel(ax,units);title(ax,titles{n});grid(ax,'on');
  if n==3,xlabel(ax,'$H$');else,xlabel(ax,'$t$');end
  draw_policy_end(ax,data{k,1}.spec,last_date,n==3);
  yline(ax,0,':','Color',[.6,.6,.6],'HandleVisibility','off');
  stabilize_neutral_axis(ax,data{k,1}.(keys{n}),data{k,2}.(keys{n}));
  if n==1,legend(ax,'Location','best','FontSize',8);end
 end
 ax=nexttile(layout);components=[data{k,1}.components;data{k,2}.components].';
 bars=bar(ax,components);for j=1:2,bars(j).FaceColor=colors(j,:);end
 xticklabels(ax,{'Hotelling','Growth','Structural'});ylabel(ax,'$100\Delta\log(\cdot)$');
 title(ax,sprintf('Extraction tilt: dates %d and %d',data{k,1}.dates));grid(ax,'on');
 yline(ax,0,':','HandleVisibility','off');stabilize_neutral_axis(ax,components(:),components(:));
 title(layout,policy_label(ids{k}));
 stems={'Figure11','Figure12'};
 save_plot(f,destination,stems{find(strcmp({'service_temporary','consumption_temporary'},ids{k}),1)},options);clear cleanup

 end
end
metric_table=cell2table(metrics,'VariableNames',{'Policy','Scenario','InitialExtractionPercent','InitialIntensityPercent', ...
 'Cumulative10Percent','Cumulative50Percent','Cumulative100Percent','CumulativeTPercent', ...
 'TotalTiltLogPoints','HotellingLogPoints','GrowthLogPoints','StructuralLogPoints', ...
 'DecompositionDate0','DecompositionDate1','StockFlowIdentityError','MaximumReturnEffectPercent'});

[baseline_table,baseline_summary]=baseline_exports(models,names,T,destination,options);

path_table=cell2table(path_rows,'VariableNames',{'Policy','Scenario','Date','ExtractionPercent','IntensityPercent', ...
 'CapitalPercent','InvestmentPercent','InvestmentRevenueSharePP','ConditionalServiceLabourSharePP', ...
 'GrossReturnPercent','CumulativeThroughDatePercent','BetaBar'});

summary=struct('source_report',report_file,'output_directory',destination,'T',T,'last_date',last_date, ...
 'policies',{ids},'metrics',metric_table,'overall_passed',report.passed,'baseline_summary',baseline_summary,'baseline_paths',baseline_table);

end

function file=find_report(source)
if isfile(source),file=source;return;end
assert(isfolder(source),'CFPlot:Source','Result source missing: %s',source);
if isfile(fullfile(source,'counterfactual_report.mat')),file=fullfile(source,'counterfactual_report.mat');return;end
files=dir(fullfile(source,'**','counterfactual_report.mat'));accepted=[];
for k=1:numel(files)
 q=load(fullfile(files(k).folder,files(k).name),'report');
 if isfield(q,'report')&&isfield(q.report,'passed')&&q.report.passed,accepted(end+1)=k;end 
end
assert(~isempty(accepted),'CFPlot:Source','No completed accepted counterfactual report found');
[~,j]=max([files(accepted).datenum]);k=accepted(j);file=fullfile(files(k).folder,files(k).name);

end
function r=get_case(model,id,T)
k=find(cellfun(@(q)strcmp(q.id,id),model.policy_cases));assert(numel(k)==1,'CFPlot:Case','Policy missing/duplicated: %s',id);
j=find(model.horizons==T);assert(numel(j)==1,'CFPlot:Horizon','Requested horizon missing');
r=model.policy_cases{k}.results{j};assert(r.passed,'CFPlot:Case','Policy equilibrium not validated');
end
function b=get_baseline(model,T)
j=find(model.horizons==T);assert(numel(j)==1,'CFPlot:Horizon','Baseline horizon missing');
b=model.baseline_cases{j};assert(b.passed,'CFPlot:Baseline','Baseline equilibrium not validated');
end
function u=policy_data(r,b,T)
p=r.path;bp=b.path;
Y=p.p_m(:).*p.c_m(:)+p.p_s(:).*p.c_s(:)+p.I(:);
Yb=bp.p_m(:).*bp.c_m(:)+bp.p_s(:).*bp.c_s(:)+bp.I(:);
u.extraction=100*(p.R_input(:)./bp.R_input(:)-1);
u.intensity=100*((p.R_input(:)./Y)./(bp.R_input(:)./Yb)-1);
u.cumulative=100*(cumsum(p.R_input(1:T))./cumsum(bp.R_input(1:T))-1);u.cumulative=u.cumulative(:);
u.K=100*(p.K(:)./bp.K(:)-1);u.I=100*(p.I(:)./bp.I(:)-1);
u.investment_share_pp=100*(p.I(:)./Y-bp.I(:)./Yb);
u.service_labour_pp=100*(p.l_s(:)./(p.l_m(:)+p.l_s(:))-bp.l_s(:)./(bp.l_m(:)+bp.l_s(:)));
q=p.investment_user_price(:);qb=bp.investment_user_price(:);delta=r.p.delta;
rnext=p.r(:);rnext=rnext(2:end);rbnext=bp.r(:);rbnext=rbnext(2:end);
ret=(rnext+(1-delta)*q(2:end))./q(1:end-1);
retb=(rbnext+(1-delta)*qb(2:end))./qb(1:end-1);
u.gross_return=100*(ret./retb-1);
u.beta_bar=(r.p.beta_m*p.p_m(:).*p.c_m(:)+r.p.beta_s*p.p_s(:).*p.c_s(:)+r.p.beta_x*p.I(:))./Y;
u.spec=r.policy_specification;
u.components=r.decomposition.component_log_points(:).';u.dates=r.decomposition.dates;
flow=sum(p.R_input(1:T))-sum(bp.R_input(1:T));stock=bp.R(T+1)-p.R(T+1);
u.initial_stock_flow_error=abs(flow-stock)/max(1,abs(stock));
assert(u.initial_stock_flow_error<1e-5,'CFPlot:StockFlow','Resource stock-flow identity failed');
end
function draw_policy_end(ax,spec,last,is_cumulative)
if spec.permanent,return;end
boundary=spec.start+spec.duration; 
if ~is_cumulative,boundary=boundary-.5;end
if boundary<=last,xline(ax,boundary,':','Color',[.6,.6,.6],'HandleVisibility','off');end
end
function stabilize_neutral_axis(ax,a,b)
if max(abs([a(:);b(:)]))<1e-6,ylim(ax,[-1e-6,1e-6]);end
end
function label=policy_label(id)
switch id
 case 'service_temporary',label='Temporary services subsidy';
 case 'service_permanent',label='Permanent services subsidy';
 case 'consumption_temporary',label='Temporary consumption tax';
 case 'consumption_permanent',label='Permanent consumption tax (neutrality check)';
 case 'carbon_temporary',label='Temporary carbon tax';
 case 'investment_purchase_20',label='Investment purchase subsidy (20%)';
 case 'investment_purchase_matched',label='Matched investment purchase subsidy (1/6)';
 otherwise,label=strrep(id,'_',' ');
end
end

function [paths,summary]=baseline_exports(models,names,T,out,o)
rows=cell(0,9);stats=cell(0,7);
for j=1:numel(models)
 b=get_baseline(models{j},T);p=b.path;
 Y=p.p_m(:).*p.c_m(:)+p.p_s(:).*p.c_s(:)+p.I(:);
 intensity=p.R_input(:)./Y;
 bar=(b.p.beta_m*p.p_m(:).*p.c_m(:)+b.p.beta_s*p.p_s(:).*p.c_s(:)+b.p.beta_x*p.I(:))./Y;
 for t=0:T
  cum=NaN;if t>0,cum=sum(p.R_input(1:t));end
  rows(end+1,:)={names{j},t,Y(t+1),p.R_input(t+1),intensity(t+1),bar(t+1),p.K(t+1),p.R(t+1),cum}; 
 end
 stats(end+1,:)={names{j},intensity(1),intensity(11),intensity(T+1),bar(1),bar(T+1),sum(p.R_input(1:T))}; 
end
paths=cell2table(rows,'VariableNames',{'Economy','Date','FinalOutput','Extraction','ResourceIntensity','BetaBar','Capital','Reserves','CumulativeHFlows'});
summary=cell2table(stats,'VariableNames',{'Economy','InitialIntensity','Date10Intensity','TerminalIntensity','InitialBetaBar','TerminalBetaBar','CumulativeTFlows'});

if ~o.make_plots,return;end
f=paper_figure(2,2);c=onCleanup(@()close(f)); 
tiledlayout(2,2);keys={'ResourceIntensity','BetaBar','Extraction','CumulativeHFlows'};
for k=1:4
 nexttile;hold on;
 for j=1:2
  d=paths(strcmp(paths.Economy,names{j}),:);plot(d.Date,d.(keys{k}),'LineWidth',o.line_width,'DisplayName',names{j});
 end
 xlabel('$t$');if k==4,xlabel('$H$');end;ylabel(math_axis(keys{k}));grid on;legend('Location','best');
end
save_plot(f,out,'Figure13',o);
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
function s=math_axis(key)
switch key
 case 'extraction',s='$100(e_t/e_t^b-1)$';
 case 'intensity',s='$100[(e_t/Y_t)/(e_t^b/Y_t^b)-1]$';
 case 'cumulative',s='$100(\mathcal{E}_H/\mathcal{E}_H^b-1)$';
 case 'K',s='$100(K_t/K_t^b-1)$';
 case 'I',s='$100(I_t/I_t^b-1)$';
 case 'investment_share_pp',s='$100\Delta(I_t/Y_t)$';
 case 'service_labour_pp',s='$100\Delta[L_{s,t}/(L_{m,t}+L_{s,t})]$';
 case 'gross_return',s='$100(\mathcal{R}_{t+1}/\mathcal{R}_{t+1}^b-1)$';
 case 'ResourceIntensity',s='$e_t/Y_t$';
 case 'BetaBar',s='$\bar\beta_t$';
 case 'Extraction',s='$e_t$';
 case 'CumulativeHFlows',s='$\mathcal{E}_H$';
 otherwise,error('Writing:Label','Unknown metric %s',key);
end
end
