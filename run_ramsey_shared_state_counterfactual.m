function report=run_ramsey_shared_state_counterfactual(w)
root=setup_solver();addpath(fullfile(root,'ramsey'));addpath(fullfile(root,'ramsey_shared'));
if nargin<1||isempty(w),w=ramsey_shared_state_config();end
if isempty(w.reference_report_file)
 z=load(fullfile(w.reference_output_root,'latest_ramsey.mat'),'entry');
 assert(z.entry.passed,'RamseyCF:Reference','Original Ramsey run must pass first');
 w.reference_report_file=z.entry.report_file;
end
z=load(w.reference_report_file,'report');reference=z.report;
assert(reference.passed&&isequal(reference.configuration.tools,{'resource'}), ...
 'RamseyCF:Reference','Requires a passed RESOURCE-ONLY original Ramsey report');
if isempty(w.climate_report_file),w.climate_report_file=reference.climate_source;end
assert(isfile(w.climate_report_file),'RamseyCF:Source','Set climate_report_file to the original accepted climate_report.mat');
z=load(w.climate_report_file,'report');original=z.report;
assert(original.passed,'RamseyCF:Source','Original climate parameterization did not pass');
assert(isequal(w.ramsey.tools,{'resource'}),'RamseyCF:Tools','This comparison is resource-only');
assert(w.ramsey.eta==reference.scenario.eta&&abs(w.ramsey.loss_target-reference.scenario.loss_target)<1e-12&& ...
 isequal(w.ramsey.resource_tax_bounds,reference.configuration.resource_tax_bounds)&& ...
 strcmp(w.ramsey.terminal_policy,reference.configuration.terminal_policy), ...
 'RamseyCF:Comparison','Match original chi scenario, resource-tax bounds and terminal continuation');
ix=find([original.main.scenarios.eta]==reference.scenario.eta& ...
 abs([original.main.scenarios.loss_target]-reference.scenario.loss_target)<1e-12,1);
assert(~isempty(ix)&&abs(original.main.scenarios(ix).chi/reference.scenario.chi-1)<1e-12, ...
 'RamseyCF:Comparison','Climate source and original Ramsey chi differ');
assert(isequal(original.main.inherited,reference.common_inherited)&& ...
 isequal(original.main.anchors,reference.common_normalization), ...
 'RamseyCF:Comparison','Original history or climate normalization differs');
m=original.configuration;base=m.base;base.make_plots=false;base.audit.derivatives=false;base.audit.restart=false;
horizons=m.T;if w.ramsey.check_horizon,horizons=[m.T m.check_T];end
assert(all(ismember(horizons,cellfun(@(a)a.calendar_terminal,reference.runs))), ...
 'RamseyCF:Horizon','Original Ramsey report missing requested horizon');
V0=[original.main.baseline.path.p_m(1)*original.main.baseline.path.c_m(1), ...
 original.main.baseline.path.p_s(1)*original.main.baseline.path.c_s(1),original.main.baseline.path.I(1)];
weights=V0/sum(V0);beta0=weights*[base.p.beta_m;base.p.beta_s;base.p.beta_x];
beta=w.common_beta;if isempty(beta),beta=beta0;end
cfg=base;cfg.horizons=horizons;cfg.output_root=w.output_root;cfg.run_kind='experiments';cfg.output_tag='ramsey_shared_state_flat_intensity';
[out,started]=create_run_directory(cfg);
report=struct('passed',false,'configuration',w,'output_dir',out,'started_at',started, ...
 'reference_file',w.reference_report_file,'original_climate_file',w.climate_report_file, ...
 'common_beta',beta,'original_initial_beta_bar',beta0,'original_initial_revenue_weights',weights, ...
 'own_history',false,'common_inherited_states',true,'production_switch_unexpected',true,'recalibrated',false,'chi_recalibrated',false,'isolated_causal_channel',false, ...
 'global_optimum_certified',false,'reference',reference);
entry=struct('passed',false,'report_file',fullfile(out,'ramsey_shared_state_report.mat'));
save(fullfile(w.output_root,'latest_ramsey_shared_state.mat'),'entry');
previous=get(0,'Diary');previous_file=get(0,'DiaryFile');diary(fullfile(out,'console.txt'));diary on;
cleanup=onCleanup(@()restore_diary(previous,previous_file)); 
try

 flatcfg=ramsey_shared_production_config(base,beta,1);
 report.flat_configuration=flatcfg;report.climate_kernel_weights=m.a;
 report.flat_baselines=cell(1,numel(horizons));report.production_histories=cell(size(horizons));
 for hi=1:numel(horizons)
  T=horizons(hi);ref=original.main.baseline;if hi>1,ref=original.check.baseline;end
  [b,h]=ramsey_shared_baseline(base,beta,ref,original.main.inherited,w,out,T);
  assert(abs(b.validation.diagnostics.terminal.investment_growth_gap)<m.baseline_tail_gap_tolerance, ...
   'RamseyCF:Tail','Flattened baseline terminal gap too large');
  report.flat_baselines{hi}=b;report.production_histories{hi}=h;
  save(fullfile(out,'ramsey_shared_state_report.mat'),'report');
 end
 b=report.flat_baselines{1};inherited=original.main.inherited;
 assert(isequal(inherited,reference.common_inherited),'RamseyShared:States','Inherited states differ');
 report.shared_state_audit=struct('passed',true,'common_K',inherited.K,'common_R',inherited.R, ...
  'common_environment_components',inherited.environment_components,'history_shared_through',m.L-1, ...
  'production_switch_date',m.L,'switch_unexpected',true);
 mc=m;mc.base=flatcfg;
 main=struct('baseline',b,'inherited',inherited,'anchors',original.main.anchors, ...
  'scenarios',original.main.scenarios,'benchmark',original.main.benchmark);
 source=struct('passed',true,'kind','shared_state_production_switch_handoff', ...
  'configuration',mc,'main',main,'check',struct(), ...
  'original_parameterization_file',w.climate_report_file,'chi_recalibrated',false, ...
  'own_history',false,'common_inherited_states',true,'production_switch_unexpected',true,'normalization_changed',false,'climate_fit_performed',false);
 if numel(horizons)>1,source.check.baseline=report.flat_baselines{2};end
 handoff=fullfile(out,'counterfactual_environment_handoff.mat');
 save_input_handoff(handoff,source);
 report.counterfactual_inherited=inherited;report.counterfactual_environment_handoff=handoff;
 o=w.ramsey;o.climate_report_file=handoff;o.output_root=fullfile(out,'counterfactual_ramsey');

 report.counterfactual=run_ramsey(o);
 report.comparison=ramsey_shared_export(reference,report.counterfactual,report,out,w);
 report.passed=report.counterfactual.passed;
 report.completed_at=char(datetime('now'));save(fullfile(out,'ramsey_shared_state_report.mat'),'report');
 entry.passed=report.passed;save(fullfile(w.output_root,'latest_ramsey_shared_state.mat'),'entry');

catch ME
 report.passed=false;report.error_identifier=ME.identifier;report.error_message=ME.message;
 save(fullfile(out,'ramsey_shared_state_report.mat'),'report');entry.passed=false;
 save(fullfile(w.output_root,'latest_ramsey_shared_state.mat'),'entry');
 rethrow(ME);
end
end
function restore_diary(state,file)
diary off;if strcmpi(state,'on'),diary(file);diary on;end
end

function save_input_handoff(file,source)
report=source;save(file,'report');
end
