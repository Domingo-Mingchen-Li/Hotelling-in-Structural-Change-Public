function report = run_counterfactual(m)
setup_solver();if nargin<1||isempty(m),m=counterfactual_config();end
assert(strcmp(m.base.model,'paper'),'CF:Mode','Counterfactual requires paper mode');
assert(m.T>=50&&m.T==fix(m.T)&&isscalar(m.check_T)&&m.check_T>m.T&&m.check_T==fix(m.check_T), ...
 'CF:Horizon','Specify T and a larger check_T');
assert(strcmp(m.share_rule,'baseline_initial_output_weighted')&& ...
 strcmp(m.adjustment,'preserve_capital_labour_ratio'),'CF:Definition','Unsupported experiment definition');
assert(m.production_initial_step>0&&m.production_initial_step<=m.production_maximum_step&& ...
 m.production_maximum_step<=1&&m.production_minimum_step>0&& ...
 m.production_minimum_step<=m.production_initial_step&&m.step_growth>1&& ...
 m.policy_initial_step>0&&m.policy_initial_step<=1&&m.policy_minimum_step>0&& ...
 m.policy_minimum_step<=m.policy_initial_step,'CF:Options','Invalid continuation options');
assert(m.maximum_continuation_attempts>=1&&m.maximum_continuation_attempts==fix(m.maximum_continuation_attempts), ...
 'CF:Options','Invalid attempt limit');
assert(all(m.cumulative_periods>0)&all(m.cumulative_periods==fix(m.cumulative_periods))& ...
 all(m.cumulative_periods<=m.T),'CF:Windows','Invalid H-period cumulative windows');
cfg=m.base;cfg.horizons=m.T;cfg.warm_start_file=m.seed_file;
assert(numel(unique(m.policy_ids))==numel(m.policy_ids)&&all(ismember(m.policy_ids,{cfg.policies.id})), ...
 'CF:Policy','Missing or duplicate policies');
cfg.policies=cfg.policies(ismember({cfg.policies.id},m.policy_ids));cfg.policy_ids=m.policy_ids;
cfg.audit.derivatives=m.audit_derivatives;cfg.audit.restart=m.audit_restart;cfg.audit.horizons=false;
cfg.newton_options.display=m.display_inner;cfg.make_plots=false;
cfg.output_root=m.output_root;cfg.output_tag='flat_resource_input_shares';cfg.run_kind='experiments';cfg=prepare_config(cfg);
[out,started]=create_run_directory(cfg);
report=struct('passed',false,'configuration',m,'output_dir',out,'started_at',started, ...
 'recalibrated',false,'physical_emissions_calibrated',false,'isolated_causal_channel',false);
save(fullfile(out,'counterfactual_configuration.mat'),'m','cfg');
entry=struct('passed',false,'report_file',fullfile(out,'counterfactual_report.mat'));
save(fullfile(m.output_root,'latest_counterfactual.mat'),'entry');
old_diary=get(0,'Diary');old_diary_file=get(0,'DiaryFile');
diary(fullfile(out,'console.txt'));diary on;
log_cleanup=onCleanup(@()restore_diary(old_diary,old_diary_file)); 

try
 assert(isfile(m.seed_file),'CF:Seed','Warm-start file missing: %s',m.seed_file);
 raw=load(m.seed_file,'warm_starts');pool=raw.warm_starts(:).';

 report.reference=counterfactual_run_cases(cfg,pool,m,out,'reference');
 b=report.reference.baseline_cases{1};V=[b.path.p_m(1)*b.path.c_m(1), ...
  b.path.p_s(1)*b.path.c_s(1),b.path.I(1)];weights=V/sum(V);
 original_betas=[cfg.p.beta_m,cfg.p.beta_s,cfg.p.beta_x];beta_initial=original_betas*weights.';
 common_beta=beta_initial;manual_override=false;
 if ~isempty(m.common_beta),common_beta=m.common_beta;manual_override=true;end
 assert(isscalar(common_beta)&&isfinite(common_beta)&&common_beta>0&&common_beta<1,'CF:Shares','Invalid common beta');
 report.common_beta=common_beta;report.reference_initial_beta_bar=beta_initial;
 report.reference_initial_revenue_weights=weights;report.manual_beta_override=manual_override;

 flat=counterfactual_production_config(cfg,common_beta,1);flat.horizons=[m.T,m.check_T];flat=prepare_config(flat);
 for sector={'m','s','x','e'}
  i=sector{1};a=flat.p.(['alpha_' i]);be=flat.p.(['beta_' i]);

 end
 g=transition_abgp_rates(flat.p,'paper');

 report.flat_configuration=flat;save(fullfile(out,'counterfactual_report.mat'),'report');

 [anchor,report.production_continuation]=counterfactual_prepare_anchor(cfg,common_beta,b,m,out);

 report.flat=counterfactual_run_cases(flat,{anchor},m,out,'flat');

 [report.channel_checks,report.channel_closed]=counterfactual_channel_audit(report.flat,common_beta,m);
 report.horizon_screen_passed=report.flat.horizon_screen_passed;
 report.cumulative_horizon_checks=counterfactual_cumulative_audit(report.flat,m);
 report.cumulative_horizon_passed=all([report.cumulative_horizon_checks.passed]);
 report.finite_equilibrium_passed=report.reference.passed&&report.flat.passed;
 report.passed=report.finite_equilibrium_passed&&report.channel_closed&& ...
  report.horizon_screen_passed&&report.cumulative_horizon_passed;
 report.completed_at=char(datetime('now'));
 save(fullfile(out,'counterfactual_report.mat'),'report');
 entry.passed=report.passed;save(fullfile(m.output_root,'latest_counterfactual.mat'),'entry');

 report.exports=counterfactual_export(report,out,m);save(fullfile(out,'counterfactual_report.mat'),'report');

 assert(report.passed,'CF:Acceptance','Counterfactual acceptance incomplete; inspect saved reports and horizon tables');
catch ME
 report.passed=false;report.error_identifier=ME.identifier;report.error_message=ME.message;
 save(fullfile(out,'counterfactual_report.mat'),'report');entry.passed=false;
 save(fullfile(m.output_root,'latest_counterfactual.mat'),'entry');

 rethrow(ME);
end
end

function restore_diary(state,file)
diary off;
if strcmpi(state,'on'),diary(file);diary on;end
end
