function report=run_industrial_welfare(w)
setup_solver();if nargin<1||isempty(w),w=industrial_welfare_config();end
if isempty(w.climate_report_file)
 marker=fullfile(w.climate_output_root,'latest_climate.mat');
 assert(isfile(marker),'Welfare:Source','Run run_climate_parameterization first, or set w.climate_report_file');
 z=load(marker,'entry');assert(z.entry.passed,'Welfare:Source','Latest climate run did not pass');
 w.climate_report_file=z.entry.report_file;
end
assert(isfile(w.climate_report_file),'Welfare:Source','Missing climate_report.mat: %s',w.climate_report_file);
z=load(w.climate_report_file,'report');source=z.report;
assert(source.passed,'Welfare:Source','Climate parameterization was not accepted');
m=source.configuration;cfg=m.base;assert(strcmp(cfg.model,'paper'),'Welfare:Mode','Requires paper mode');
assert(~isempty(w.policies)&&numel(unique({w.policies.id}))==numel(w.policies),'Welfare:Policies','Empty or duplicate policy IDs');
for j=1:numel(w.policies)
 s=w.policies(j);
 assert(s.start==0&&s.duration>=2&&s.duration==fix(s.duration)&&s.duration<m.T-m.L, ...
  'Welfare:Timing','Policies must start at announcement and last >=2 nodes');
 assert(s.tau_e==0&&s.tau_int==0,'Welfare:Scope','This exercise excludes resource/intermediate taxes');
 assert(~(s.permanent&&s.tau_x~=0),'Welfare:Tail','Permanent investment subsidies need a separate closure');
end
assert(w.horizon_anchor_tolerance>0&&w.sign_anchor_tolerance>0&&w.neutrality_tolerance>0, ...
 'Welfare:Tolerance','Require positive tolerances');
cfg.horizons=m.T;if w.check_horizon,cfg.horizons=[m.T m.check_T];end
cfg.output_root=w.output_root;cfg.run_kind='experiments';cfg.output_tag='industrial_welfare';
cfg.make_plots=false;cfg.newton_options.display=w.display_inner;
cfg.audit.derivatives=w.audit_derivatives;cfg.audit.restart=w.audit_restart;
[out,started]=create_run_directory(cfg);
report=struct('passed',false,'configuration',w,'climate_source',w.climate_report_file, ...
 'started_at',started,'output_dir',out,'policy_optimized',false,'physical_emissions_calibrated',false, ...
 'announcement_unexpected',true,'source_parameters',source.main.scenarios, ...
 'common_inherited',source.main.inherited,'common_normalization',source.main.anchors);
entry=struct('passed',false,'report_file',fullfile(out,'industrial_welfare_report.mat'));
save(fullfile(w.output_root,'latest_industrial_welfare.mat'),'entry');
previous=get(0,'Diary');previous_file=get(0,'DiaryFile');diary(fullfile(out,'console.txt'));diary on;
cleanup=onCleanup(@()restore_diary(previous,previous_file)); 
try

 horizons=m.T;if w.check_horizon,horizons=[m.T m.check_T];end
 rows=[];report.runs=cell(1,numel(horizons));
 local=m;local.announcement_calendar_date=m.L;local.L=0;
 local.initial_components=source.main.inherited.environment_components;
 anchor=source.main.benchmark.utility_loss;
 for hi=1:numel(horizons)
  T=horizons(hi);N=T-m.L;
  if hi==1,old=source.main.baseline;else,old=source.check.baseline;end
  warm=shift_seed(old,m.L,source.main.inherited);
  cc=cfg;cc.targets.K0=source.main.inherited.K;cc.targets.R0=source.main.inherited.R;
  cc.A0=source.main.inherited.A;cc.scale_h0=warm.path.h(1);cc.scale_E0=warm.path.E(1);
  cc=prepare_config(cc);

  zero=w.policies(1);zero.id='baseline';zero.permanent=false;
  for f={'tau_c_m','tau_c_s','tau_e','tau_int','tau_x'},zero.(f{1})=0;end
  zero.investment_policy_kind='none';
  b=solve_transition(cc,zero,N,warm,[]);
  assert(b.passed,'Welfare:Baseline','Restarted baseline failed');
  assert(abs(b.validation.diagnostics.terminal.investment_growth_gap)<m.baseline_tail_gap_tolerance, ...
   'Welfare:Tail','Restarted baseline not sufficiently close to ABGP');
  run=struct('calendar_terminal',T,'baseline',b,'cases',{cell(1,numel(w.policies))},'rows',[]);
  for j=1:numel(w.policies)
   spec=w.policies(j);
   r=solve_transition(cc,spec,N,b,b);
   assert(r.passed,'Welfare:Policy','Policy failed: %s',spec.id);
   assert(abs(r.validation.diagnostics.terminal.investment_growth_gap)<m.baseline_tail_gap_tolerance, ...
    'Welfare:Tail','Policy %s terminal gap too large',spec.id);
   a=industrial_welfare_components(r,b,local,source.main.anchors,source.main.scenarios,w);
   ne=neutrality(r,b,a,anchor,w);
   run.cases{j}=struct('solution',r,'welfare',a,'neutrality_audit',ne);
   if isempty(run.rows),run.rows=a.rows;else,run.rows=[run.rows a.rows];end 
   report.runs{hi}=run;save(fullfile(out,'industrial_welfare_report.mat'),'report');
   export_path(r,b,a,local,out,T);
  end
  report.runs{hi}=run;
  if isempty(rows),rows=run.rows;else,rows=[rows run.rows];end 
  writetable(struct2table(rows),fullfile(out,'welfare_all_horizons.csv'));
 end
 report.horizon_audit=[];
 if w.check_horizon
  report.horizon_audit=industrial_welfare_horizon(report.runs{1}.rows,report.runs{2}.rows,anchor,w);
  writetable(struct2table(report.horizon_audit),fullfile(out,'welfare_horizon_audit.csv'));
 end
 main=report.runs{1}.rows;
 writetable(struct2table(main),fullfile(out,'welfare_all_scenarios.csv'));
 ids=[main.eta]==m.benchmark_eta&[main.loss_target]==m.benchmark_loss;
 writetable(struct2table(main(ids)),fullfile(out,'welfare_benchmark.csv'));

 report.passed=~w.check_horizon||all([report.horizon_audit.passed]);
 report.completed_at=char(datetime('now'));save(fullfile(out,'industrial_welfare_report.mat'),'report');
 entry.passed=report.passed;save(fullfile(w.output_root,'latest_industrial_welfare.mat'),'entry');

 assert(report.passed,'Welfare:Horizon','Horizon/sign check failed; inspect saved CSVs before interpreting welfare');
catch ME
 report.passed=false;report.error_identifier=ME.identifier;report.error_message=ME.message;
 save(fullfile(out,'industrial_welfare_report.mat'),'report');entry.passed=false;
 save(fullfile(w.output_root,'latest_industrial_welfare.mat'),'entry');
 rethrow(ME);
end
end
function warm=shift_seed(old,L,inherited)
warm=old;warm.case_name='shifted_seed';warm.T=old.T-L;
for f=fieldnames(old.path).'
 v=old.path.(f{1});assert(isvector(v)&&numel(v)==old.T+1,'Welfare:Seed','Unexpected path shape');
 warm.path.(f{1})=v(L+1:end);
end
warm.Targets.K0=inherited.K;warm.Targets.R0=inherited.R;
end
function a=neutrality(r,b,welfare,anchor,w)
s=r.policy_specification;
applicable=s.permanent&&s.tau_c_m==s.tau_c_s&&s.tau_x==0;
a=struct('applicable',applicable,'passed',true,'allocation_error',NaN,'welfare_error_in_anchor_units',NaN);
if ~applicable,return;end
err=0;
for f={'K','R','r','h','I','R_input','c_m','c_s'},err=max(err,max(abs(r.path.(f{1})./b.path.(f{1})-1)));end
we=max(abs([welfare.rows.consumption_utility_change welfare.rows.environment_welfare_change]))/anchor;
a.allocation_error=err;a.welfare_error_in_anchor_units=we;a.passed=err<w.neutrality_tolerance&&we<w.neutrality_tolerance;

assert(a.passed,'Welfare:Neutrality','Permanent uniform consumption tax neutrality failed');
end
function export_path(r,b,a,m,out,T)
N=r.T;local_date=(0:N)';calendar_date=local_date+m.announcement_calendar_date;
weight=r.p.beta.^local_date;du=a.private_flow_difference(:);
ix=find([a.rows.eta]==m.benchmark_eta&[a.rows.loss_target]==m.benchmark_loss,1);de=a.environment_flow_difference(ix,:)';
included=double(local_date<N);cu=cumsum(weight.*du.*included);ce=cumsum(weight.*de.*included);
t=table(calendar_date,local_date,r.path.R_input(:),b.path.R_input(:),a.policy_environment.M(:), ...
 a.baseline_environment.M(:),du,de,cu,ce,cu+ce, ...
 'VariableNames',{'calendar_date','local_date','policy_extraction','baseline_extraction','policy_burden', ...
 'baseline_burden','consumption_utility_flow_difference','environment_welfare_flow_difference', ...
 'finite_consumption_cumulative','finite_environment_cumulative','finite_net_cumulative'});
writetable(t,fullfile(out,sprintf('path_%s_T%d.csv',r.case_name,T)));
end
function restore_diary(state,file)
diary off;if strcmpi(state,'on'),diary(file);diary on;end
end
