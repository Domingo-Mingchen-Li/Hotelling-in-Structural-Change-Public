function table2=export_paper_cumulative_table(source,out)
a=load(fullfile(source,'report.mat'),'report');r=a.report;
k=find(r.horizons==200,1);b=r.baseline_cases{k}.path;
H=[10 50 100 200];names=cell(numel(r.policy_cases),1);v=zeros(numel(names),numel(H));
for j=1:numel(names)
 q=r.policy_cases{j};names{j}=q.id;p=q.results{k}.path;
 for n=1:numel(H),ix=1:H(n);v(j,n)=100*(sum(p.R_input(ix))/sum(b.R_input(ix))-1);end
end
table2=array2table(v,'VariableNames',{'H10','H50','H100','H200'});
table2=addvars(table2,string(names),'Before',1,'NewVariableNames','Policy');
writetable(table2,fullfile(out,'Table2_cumulative_extraction.csv'));
fid=fopen(fullfile(out,'Table2_cumulative_extraction.tex'),'w');assert(fid>=0);done=onCleanup(@()fclose(fid)); 
fprintf(fid,'\\begin{tabular}{lrrrr}\n\\toprule\nPolicy & $H=10$ & $H=50$ & $H=100$ & $H=200$ \\\\ \n\\midrule\n');
for j=1:numel(names)
 label=strrep(names{j},'_','\_');fprintf(fid,'%s & %.6f & %.6f & %.6f & %.6f \\\\ \n',label,v(j,:));
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n');
end
