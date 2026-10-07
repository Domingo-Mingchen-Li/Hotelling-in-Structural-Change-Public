function report=run_ramsey(w)
root=setup_solver();addpath(fullfile(root,'ramsey'));
if nargin<1||isempty(w),w=ramsey_config();end
if isempty(w.climate_report_file)
 marker=fullfile(w.climate_output_root,'latest_climate.mat');
 assert(isfile(marker),'Ramsey:Source','Run climate parameterization first or set climate_report_file');
 z=load(marker,'entry');assert(z.entry.passed,'Ramsey:Source','Latest climate run did not pass');
 w.climate_report_file=z.entry.report_file;
end
z=load(w.climate_report_file,'report');source=z.report;
assert(source.passed,'Ramsey:Source','Requires accepted climate parameters');
assert(iscell(w.tools)&&~isempty(w.tools)&&numel(unique(w.tools))==numel(w.tools)&& ...
 all(ismember(w.tools,{'resource','services'})),'Ramsey:Tools','Select resource, services or both');
assert(strcmp(w.terminal_policy,'zero_from_T'),'Ramsey:Tail','Only explicit zero-from-T continuation is implemented');
for b={w.resource_tax_bounds,w.service_wedge_bounds}
 assert(numel(b{1})==2&&all(isfinite(b{1}))&&b{1}(1)>-1&&b{1}(1)<=0&&b{1}(2)>=0&&b{1}(2)>b{1}(1), ...
  'Ramsey:Bounds','Bounds must contain zero and have positive gross wedges');
end
assert(w.max_iterations>=1&&w.memory>=1&&w.max_backtracks>=1&&w.derivative_step>0&& ...
 w.projected_gradient_tolerance>0&&~isempty(w.start_amplitudes)&&all(w.start_amplitudes>=0)&& ...
 w.start_decay_years>0,'Ramsey:Options','Invalid optimization options');
m=source.configuration;sc=source.main.scenarios;
ix=find([sc.eta]==w.eta&abs([sc.loss_target]-w.loss_target)<1e-12,1);
assert(~isempty(ix),'Ramsey:Scenario','Selected chi scenario missing; re-run climate parameterization to add it');
sc=sc(ix);cfg=m.base;cfg.output_root=w.output_root;cfg.run_kind='experiments';cfg.output_tag=['ramsey_' strjoin(w.tools,'_')];
horizons=m.T;if w.check_horizon,horizons=[m.T m.check_T];end
cfg.horizons=horizons;cfg.make_plots=false;cfg.audit.derivatives=false;cfg.audit.restart=false;
[out,started]=create_run_directory(cfg);
report=struct('passed',false,'configuration',w,'climate_source',w.climate_report_file, ...
 'output_dir',out,'started_at',started,'scenario',sc,'common_inherited',source.main.inherited, ...
 'common_normalization',source.main.anchors,'runs',{cell(1,numel(horizons))}, ...
 'terminal_policy','zero from calendar T onward','global_optimum_certified',false, ...
 'domain','strictly positive scarcity rents and exhaustion-type ABGP continuation', ...
 'announcement_unexpected',true,'first_best_certified',false,'zero_rent_corners_included',false);
entry=struct('passed',false,'report_file',fullfile(out,'ramsey_report.mat'));
save(fullfile(w.output_root,'latest_ramsey.mat'),'entry');
previous=get(0,'Diary');previous_file=get(0,'DiaryFile');diary(fullfile(out,'console.txt'));diary on;
cleanup=onCleanup(@()restore_diary(previous,previous_file)); 
try

 summary=[];
 for hi=1:numel(horizons)
  T=horizons(hi);N=T-m.L;old=source.main.baseline;if hi>1,old=source.check.baseline;end
  warm=old;warm.T=N;warm.case_name='ramsey_shifted_seed';
  for f=fieldnames(old.path).',warm.path.(f{1})=old.path.(f{1})(m.L+1:end);end
  warm.Targets.K0=source.main.inherited.K;warm.Targets.R0=source.main.inherited.R;
  cc=cfg;cc.targets=warm.Targets;cc.A0=source.main.inherited.A;
  cc.scale_h0=warm.path.h(1);cc.scale_E0=warm.path.E(1);cc.newton_options=w.inner_options;
  cc.newton_options.display=w.display_inner;cc=prepare_config(cc);
  zero=struct('id','baseline','start',0,'duration',10,'permanent',false, ...
   'tau_c_m',0,'tau_c_s',0,'tau_e',0,'tau_int',0,'tau_x',0,'investment_policy_kind','none');

  b=solve_transition(cc,zero,N,warm,[]);assert(b.passed,'Ramsey:Baseline','Restart baseline failed');
  ctx=struct('problem',b.problem,'tools',{w.tools},'options',w,'scenario',sc, ...
   'anchors',source.main.anchors,'anchor',source.main.benchmark.utility_loss,'baseline_value',0);
  ctx.local_climate=m;ctx.local_climate.L=0;ctx.local_climate.initial_components=source.main.inherited.environment_components;
  ctx.lower=[];ctx.upper=[];
  for k=1:numel(w.tools)
   bounds=w.resource_tax_bounds;if strcmp(w.tools{k},'services'),bounds=w.service_wedge_bounds;end
   ctx.lower=[ctx.lower;log1p(bounds(1))*ones(N,1)];ctx.upper=[ctx.upper;log1p(bounds(2))*ones(N,1)]; 
  end
  bv=ramsey_value_gradient(b.x,b.problem,ctx,false);ctx.baseline_value=bv.value;
  e0=ramsey_evaluate(zeros(numel(ctx.lower),1),ctx,b.x,true);assert(e0.success,'Ramsey:Baseline','Initial gradient failed');
  run=struct('calendar_terminal',T,'baseline',b,'baseline_welfare',bv, ...
   'initial_gradient_audit',[],'starts',{cell(1,numel(w.start_amplitudes))},'passed',false);
  if w.audit_gradient,run.initial_gradient_audit=ramsey_gradient_audit(e0,ctx);end
  for j=1:numel(w.start_amplitudes)
   u=[];
   for k=1:numel(w.tools)
    amplitude=w.start_amplitudes(j);if strcmp(w.tools{k},'services'),amplitude=-amplitude;end
    amplitude=max(expm1(ctx.lower((k-1)*N+1)),min(expm1(ctx.upper((k-1)*N+1)),amplitude));
    u=[u;log1p(amplitude)*exp(-(0:N-1)'/w.start_decay_years)]; 
   end
   if hi>1&&j==1
    oldbest=report.runs{1}.best.evaluation.u;oldN=report.runs{1}.baseline.T;u=[];
    for k=1:numel(w.tools),u=[u;oldbest((k-1)*oldN+(1:oldN));zeros(N-oldN,1)];end 
   end

   cp=fullfile(out,sprintf('checkpoint_T%d_start%d.mat',T,j));
   o=ramsey_projected_lbfgs(u,ctx,b.x,cp);run.starts{j}=o;
   report.runs{hi}=run;save(fullfile(out,'ramsey_report.mat'),'report');
  end
  gains=cellfun(@(a)a.evaluation.gain,run.starts);[~,bestid]=max(gains);run.best=run.starts{bestid};
  run.multistart_gain_spread=max(gains)-min(gains);
  run.multistart_passed=all(cellfun(@(a)a.passed,run.starts))&&run.multistart_gain_spread<w.multistart_gain_tolerance;
  e=run.best.evaluation;v=transition_validate_solution(e.x,e.problem,cc.tolerance);
  assert(v.passed,'Ramsey:Validation','Final competitive equilibrium failed independent validation');
  run.final_validation=v;run.final_fiscal_foc_audit=audit_policy_equations(transition_export_path(e.x,e.problem),e.problem);
  if w.audit_equilibrium_derivatives,run.equilibrium_derivative_audit=audit_sparse_derivatives(e.x,e.problem);end
  if w.audit_gradient,run.final_gradient_audit=ramsey_gradient_audit(e,ctx);end
  run.upper_bound_active=any(e.u>=ctx.upper-1e-6);
  run.terminal_gap_passed=abs(e.terminal_gap)<m.baseline_tail_gap_tolerance;
  run.passed=run.best.passed&&run.multistart_passed&&run.terminal_gap_passed&&e.gain>=-1e-8;

  report.runs{hi}=run;row=ramsey_export_run(run,ctx,m.L,out,w);
  if isempty(summary),summary=row;else,summary(end+1)=row;end 
  writetable(struct2table(summary),fullfile(out,'ramsey_summary.csv'));
  save(fullfile(out,'ramsey_report.mat'),'report');
 end
 report.horizon_audit=struct('applicable',w.check_horizon,'passed',true);
 if w.check_horizon
  a=report.runs{1};b=report.runs{2};D=min(w.horizon_comparison_dates,a.baseline.T);diffs=[];
  for k=1:numel(w.tools)
   name='tau_e_path';if strcmp(w.tools{k},'services'),name='tau_c_s_path';end
   diffs=[diffs;abs(a.best.evaluation.problem.policy.(name)(1:D)-b.best.evaluation.problem.policy.(name)(1:D)).']; 
  end
  dg=abs(a.best.evaluation.gain-b.best.evaluation.gain);
  report.horizon_audit=struct('applicable',true,'gain_difference_in_anchor_units',dg, ...
   'maximum_early_wedge_difference',max(diffs),'comparison_local_dates',D, ...
   'passed',dg<w.horizon_gain_tolerance&&max(diffs)<w.horizon_early_wedge_tolerance);
  writetable(struct2table(report.horizon_audit),fullfile(out,'ramsey_horizon_audit.csv'));

 end
 report.passed=all(cellfun(@(a)a.passed,report.runs))&&report.horizon_audit.passed;
 report.completed_at=char(datetime('now'));save(fullfile(out,'ramsey_report.mat'),'report');
 entry.passed=report.passed;save(fullfile(w.output_root,'latest_ramsey.mat'),'entry');

 assert(report.passed,'Ramsey:Acceptance','Inspect checkpoints and CSVs: numerical acceptance checks failed');
catch ME
 report.passed=false;report.error_identifier=ME.identifier;report.error_message=ME.message;
 save(fullfile(out,'ramsey_report.mat'),'report');entry.passed=false;
 save(fullfile(w.output_root,'latest_ramsey.mat'),'entry');
 rethrow(ME);
end
end
function restore_diary(state,file)
diary off;if strcmpi(state,'on'),diary(file);diary on;end
end
