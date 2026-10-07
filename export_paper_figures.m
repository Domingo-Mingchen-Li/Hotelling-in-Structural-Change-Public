function output_dir = export_paper_figures(experiment_dir,robustness_dir,output_dir,mode)
root=fileparts(mfilename('fullpath'));
if nargin<4,mode='main';end
if nargin<1||isempty(experiment_dir)
 experiment_dir=fullfile(root,'results','paper','transition');
end
if nargin<2
 robustness_dir=fullfile(root,'results','paper','horizon');
end
raw=load(fullfile(experiment_dir,'report.mat'),'report');report=raw.report;
assert(report.passed,'fig:Validation','Experiment report must pass');
T=200;k=find(report.horizons==T,1);
assert(~isempty(k),'fig:Horizon','A T=200 result is required');
b=report.baseline_cases{k};assert(b.passed,'fig:Baseline','Invalid baseline');b=b.path;
last=40;idx=1:last+1;t=0:last;
if nargin<3||isempty(output_dir),output_dir=fullfile(root,'outputs','paper_transition');end
if ~isfolder(output_dir),mkdir(output_dir);end
colors=[0.08 0.28 0.48;0.74 0.27 0.16;0.20 0.48 0.31;0.47 0.32 0.60];
if strcmp(mode,'main')
Yb=output(b);Lb=labour_shares(b);
st=solution('service_temporary');sp=solution('service_permanent');
ct=solution('consumption_temporary');carbon=solution('carbon_temporary');
ip=solution('investment_purchase_20');im=solution('investment_purchase_matched');
assert(st.policy_specification.duration==10,'fig:Duration','Expected ten policy nodes');
f=newfigure(2,4);
y={b.K,b.p_m./b.p_s,100*Lb(1,:),100*Lb(2,:),b.R/1e5,b.R_input/1e4,b.R_input./Yb, ...
 100*b.l_s./(b.l_m+b.l_s)};
titles={'Capital stock','Relative price p_m/p_s','Manufacturing labour share', ...
 'Services labour share','Remaining resource stock','Extraction flow','Resource intensity e/Y', ...
 {'Services share within','manufacturing and services'}};
units={'Model units','Price ratio','Percent of total labour','Percent of total labour', ...
 'Model units / 100,000','Model units / 10,000','Model units / nominal output','Percent'};
for j=1:8,nexttile;lineplot(y{j},titles{j},units{j},colors(1,:));end
savefigure(f,'Figure1');
p=st.path;Y=output(p);L=labour_shares(p);f=newfigure(2,3);
y={effect(p.R,b.R),effect(p.I,b.I),100*(p.I./Y-b.I./Yb), ...
 100*(L(2,:)-Lb(2,:)),effect(p.R_input,b.R_input),effect(p.R_input./Y,b.R_input./Yb)};
titles={'Remaining resource stock','Investment','Investment / nominal output', ...
 'Services labour share','Extraction flow','Resource intensity e/Y'};
units={'Deviation (%)','Deviation (%)','Difference (pp)','Difference (pp)','Deviation (%)','Deviation (%)'};
for j=1:6,nexttile;lineplot(y{j},titles{j},units{j},colors(2,:));policyline(st);end
savefigure(f,'Figure2');
f=newfigure(1,1);nexttile;hold on;
plot(t,b.R_input(idx)/1e4,'--','Color',colors(1,:),'LineWidth',1.6);
plot(t,p.R_input(idx)/1e4,'-','Color',colors(2,:),'LineWidth',1.6);
format_axis('Extraction flow','Model units / 10,000');policyline(st);
legend({'Baseline','Temporary services subsidy'},'Location','northeast','Box','off');
savefigure(f,'Figure3');
f=newfigure(1,1);nexttile;decomposition(st);savefigure(f,'Figure4');
f=newfigure(1,1);nexttile;decomposition(sp);savefigure(f,'Figure5');
p=ct.path;Y=output(p);f=newfigure(2,3);nexttile;decomposition(ct);
y={effect(p.R_input,b.R_input),effect(p.R,b.R),effect(p.K,b.K), ...
 100*(p.I./Y-b.I./Yb),effect(p.I,b.I)};
titles={'Extraction flow','Remaining resource stock','Capital stock','Investment / nominal output','Investment'};
for j=1:5
 nexttile;if j==4,u='Difference (pp)';else,u='Deviation (%)';end
 lineplot(y{j},titles{j},u,colors(3,:));policyline(ct);
end
savefigure(f,'Figure6');
p=carbon.path;L=labour_shares(p);f=newfigure(2,2);
titles={'Manufacturing labour share','Services labour share','Investment labour share','Resource-processing labour share'};
for j=1:4,nexttile;lineplot(100*(L(j,:)-Lb(j,:)),titles{j},'Difference (pp)',colors(j,:));policyline(carbon);end
savefigure(f,'Figure9');
cases={ct,ip,im};labels={'Consumption tax 20%','Investment purchase subsidy 20%', ...
 'Investment purchase subsidy 1/6'};f=newfigure(1,3);
for j=1:3
 nexttile;hold on;
 for n=1:3
  p=cases{n}.path;
  if j==1,y=effect(p.R_input,b.R_input);titletext='Extraction flow';
  elseif j==2,y=effect(p.I,b.I);titletext='Investment';
  else,y=effect(p.R_input./output(p),b.R_input./Yb);titletext='Resource intensity e/Y';end
  plot(t,y(idx),'Color',colors(n,:),'LineWidth',1.6);
 end
 format_axis(titletext,'Deviation (%)');policyline(ct);
 if j==1
   lg=legend(labels,'Box','off','NumColumns',1);lg.Layout.Tile='south';
  end
end
savefigure(f,'Figure7');
cases={ct,ip,im};values=zeros(3,3);dates=ct.decomposition.dates;
for n=1:3
 d=cases{n}.decomposition;
 assert(isequal(d.dates,dates),'fig:DecompositionDates', ...
  'All compared decompositions must use the same endpoint dates');
 values(n,:)=reshape(d.component_log_points,1,3);
 assert(all(isfinite(values(n,:)))&& ...
  abs(sum(values(n,:))-d.total_log_points)<1e-7, ...
  'fig:DecompositionIdentity','Invalid accounting decomposition');
end
f=newfigure(1,1);set(f,'Position',[2,2,23,12]);nexttile;
bars=bar(1:3,values,'grouped');palette=colors([2,4,3],:);
for j=1:3,bars(j).FaceColor=palette(j,:);end
xticks(1:3);xticklabels({'Consumption tax 20%', ...
 'Investment subsidy 20%','Investment subsidy 16.67%'});
ylabel('Log points (100 x log deviation)');
title(sprintf('Extraction-ratio decomposition: t=%d versus t=%d',dates), ...
 'FontWeight','normal','Interpreter','none');
yline(0,'-','Color',[.5 .5 .5],'HandleVisibility','off');grid on;box off;
ax=gca;ax.YAxis.Exponent=0;
set(ax,'FontName','Times New Roman','FontSize',10,'TickDir','out');
span=max(max(values(:))-min(values(:)),1);
ylim([min([0;values(:)])-.15*span,max([0;values(:)])+.18*span]);
xlim([.5,3.5]);
lg=legend(bars,{'Hotelling','Growth','Structure'},'Box','off','NumColumns',3);
lg.Layout.Tile='south';drawnow;
for j=1:3
 for n=1:3
  v=values(n,j);
  if v>=0,align='bottom';offset=.02*span;else,align='top';offset=-.02*span;end
  text(bars(j).XEndPoints(n),v+offset,sprintf('%.2f',v), ...
   'HorizontalAlignment','center','VerticalAlignment',align,'FontSize',10);
 end
end
savefigure(f,'Figure8');

cases={st,sp,ct,carbon};labels={'Temporary services subsidy','Permanent services subsidy', ...
 'Temporary consumption tax','Temporary carbon tax'};f=newfigure(1,1);nexttile;hold on;
for n=1:4
 y=effect(cumsum(cases{n}.path.R_input),cumsum(b.R_input));
 plot(t,y(idx),'Color',colors(n,:),'LineWidth',1.6);
end
format_axis('Cumulative extraction from date 0','Deviation vs baseline (%)');
policyline(st);
lg=legend(labels,'Box','off','NumColumns',1);lg.Layout.Tile='south';
savefigure(f,'Figure10');
end
if strcmp(mode,'appendix')&&~isempty(robustness_dir)
 raw=load(fullfile(robustness_dir,'report.mat'),'report');hr=raw.report;
 assert(hr.passed,'fig:HorizonValidation','Robustness report must pass');
 Ts=hr.horizons;[~,ref]=max(Ts);hb=hr.baseline_cases{ref}.path;
 hc=lines(numel(Ts));hl=arrayfun(@(v)sprintf('T = %d',v),Ts,'UniformOutput',false);
 f=newfigure(1,3);vars={'K','R_input','E'};titles={'Capital','Extraction','Consumption expenditure'};
 for j=1:3
  nexttile;hold on;
  for n=1:numel(Ts)
   a=hr.baseline_cases{n}.path;av=a.(vars{j});bv=hb.(vars{j});
   y=effect(av(idx),bv(idx));
   plot(t,y,'Color',hc(n,:),'LineWidth',1.4);
  end
  format_axis(titles{j},sprintf('Deviation vs T=%d (%%)',Ts(ref)));
  if j==1,legend(hl,'Location','best','Box','off');end
 end
 savefigure(f,'FigureA1');
 f=newfigure(2,2);ids={'service_temporary','service_permanent'};
 for j=1:2
  pc=hr.policy_cases{find(cellfun(@(v)strcmp(v.id,ids{j}),hr.policy_cases),1)};
  nexttile;hold on;
  for n=1:numel(Ts)
   a=pc.results{n}.path;bb=hr.baseline_cases{n}.path;y=effect(a.R_input,bb.R_input);
   plot(t,y(idx),'Color',hc(n,:),'LineWidth',1.4);
  end
  format_axis(strrep(ids{j},'_',' '),'Extraction deviation (%)');policyline(pc.results{ref});
  if j==1,legend(hl,'Location','best','Box','off');end
  nexttile;values=zeros(numel(Ts),3);
  for n=1:numel(Ts),values(n,:)=pc.results{n}.decomposition.component_log_points;end
  bars=bar(Ts,values,'grouped');palette=colors([2,4,3],:);
  for component=1:3,bars(component).FaceColor=palette(component,:);end
  xlabel('Horizon T');ylabel('Log points');title('Extraction-ratio decomposition');grid on;box off;
  if j==1
   lg=legend({'Hotelling','Growth','Structure'},'Box','off','NumColumns',3);lg.Layout.Tile='south';
  end
 end
 savefigure(f,'FigureA2');
end

 function s=solution(id)
  j=find(cellfun(@(v)strcmp(v.id,id),report.policy_cases),1);
  assert(~isempty(j),'fig:Case','Missing policy %s',id);
  s=report.policy_cases{j}.results{k};assert(s.passed,'fig:Case','Policy did not pass');
 end
 function f=newfigure(rows,columns)
  f=figure('Visible','off','Color','w','Units','centimeters', ...
   'Position',[2,2,7*columns+2,5.5*rows+2]);
  tiledlayout(rows,columns,'TileSpacing','compact','Padding','compact');
 end
 function lineplot(y,titletext,unit,color)
  plot(t,y(idx),'Color',color,'LineWidth',1.6);format_axis(titletext,unit);
 end
 function format_axis(titletext,unit)
  title(titletext,'FontWeight','normal','Interpreter','none');xlabel('Model year t');ylabel(unit);
  grid on;box off;set(gca,'FontName','Times New Roman','FontSize',10,'TickDir','out');
  if contains(unit,'Deviation')||contains(unit,'Difference'),yline(0,':','Color',[.55 .55 .55],'HandleVisibility','off');end
  xlim([0,last]);ax=gca;ax.YAxis.Exponent=0;
 end
 function policyline(s)
  if ~s.policy_specification.permanent
   xline(s.policy_specification.start+s.policy_specification.duration-0.5,':', ...
    'Color',[.5 .5 .5],'HandleVisibility','off');
  end
 end
 function decomposition(s)
  d=s.decomposition;v=[d.total_log_points,d.component_log_points];
  h=bar(v,'FaceColor','flat');h.CData=[colors(1,:);colors(2,:);colors(4,:);colors(3,:)];
  xticks(1:4);xticklabels({'Total','Hotelling','Growth','Structure'});
  ylabel('Log points (100 x log deviation)');
  title(sprintf('Extraction-ratio tilt: t=%d versus t=%d',d.dates),'FontWeight','normal');
  yline(0,'-','Color',[.5 .5 .5]);grid on;box off;
  set(gca,'FontName','Times New Roman','FontSize',10,'TickDir','out');
  range=max(max(v)-min(v),1);ylim([min([0,v])-.15*range,max([0,v])+.18*range]);
  for q=1:4
   if v(q)>=0,align='bottom';offset=.03*range;else,align='top';offset=-.03*range;end
   text(q,v(q)+offset,sprintf('%.2f',v(q)),'HorizontalAlignment','center','VerticalAlignment',align,'FontSize',10);
  end
 end
 function savefigure(f,name)
  cleanup=onCleanup(@()close(f));drawnow;
  exportgraphics(f,fullfile(output_dir,[name '.pdf']),'ContentType','vector','BackgroundColor','white');
  exportgraphics(f,fullfile(output_dir,[name '.png']),'Resolution',300,'BackgroundColor','white');

 end
end
function y=output(p)
y=p.p_m.*p.c_m+p.p_s.*p.c_s+p.I;
assert(all(isfinite(y))&&all(y>0),'fig:Output','Invalid nominal final output');
end
function L=labour_shares(p)
L=[p.l_m(:).';p.l_s(:).';p.l_x(:).';p.l_e(:).'];
assert(all(isfinite(L(:)))&&all(L(:)>=0)&&all(sum(L,1)>0),'fig:Labour','Invalid labour inputs');
L=L./sum(L,1);
end
function y=effect(a,b)
assert(isvector(a)&&isvector(b)&&numel(a)==numel(b), ...
 'fig:EffectLength','Compare equal-length vectors at matching calendar dates');
a=a(:).';b=b(:).';
assert(all(isfinite(a))&&all(isfinite(b))&&all(b>0),'fig:Effect','Invalid comparison');
y=100*(a./b-1);
end
