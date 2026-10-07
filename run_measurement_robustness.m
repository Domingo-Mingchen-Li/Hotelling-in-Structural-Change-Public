function reports = run_measurement_robustness(case_id,m)
setup_solver();
if nargin<2||isempty(m),m=measurement_robustness_config();end
if nargin<1||isempty(case_id),ids=m.case_ids;else,ids={char(case_id)};end
configs=cell(size(ids));entries=cell(size(ids));fits=cell(size(ids));
for i=1:numel(ids),[configs{i},entries{i},fits{i}]=measurement_load_fit(m,ids{i});end
reference=m.base;reference.horizons=m.T;reference.policy_ids=m.policy_ids;
assert(all(ismember(m.policy_ids,{reference.policies.id})),'measurement:Policy','Unknown active policy');
reference.policies=reference.policies(ismember({reference.policies.id},m.policy_ids));
reference.warm_start_file=m.seed_file;reference.output_root=fullfile(m.output_root,'reference');
reference.run_kind='experiments';reference.output_tag='measurement_reference';
reference.make_plots=m.make_plots;reference.audit.derivatives=m.audit_derivatives;
reference.audit.restart=m.audit_restart;reference.audit.horizons=false;
ref=run_solver(reference);
mv=measurement_moments(ref.baseline_cases{1}.path);
assert(all(abs(mv-m.calibration.targets)<=m.calibration.fit_tolerances), ...
 'measurement:ReferenceFit','Paper reference no longer matches the configured moments; recalibrate reference first');
reports=cell(size(ids));
for i=1:numel(ids)
 cfg=configs{i};cfg.output_root=fullfile(m.output_root,ids{i},'policies');
 cfg.run_kind='experiments';cfg.output_tag=['measurement_' ids{i}];
 alt=run_solver(cfg);
 errors=measurement_moments(alt.baseline_cases{1}.path)-fits{i}.configuration.targets;
 assert(all(abs(errors)<=fits{i}.configuration.fit_tolerances), ...
  'measurement:FitDrift','Re-solved alternative baseline no longer matches calibration');
 comparison=measurement_export_comparison(ref,alt,m,ids{i});
 report=struct('passed',ref.passed&&alt.passed,'case_id',ids{i}, ...
  'configuration',m,'calibration_entry',entries{i},'calibration_report',fits{i}, ...
  'reference_output_dir',ref.output_dir,'alternative_output_dir',alt.output_dir, ...
  'comparison',comparison,'alternative_moment_errors',errors, ...
  'economic_sign_robustness_certified',false,'infinite_horizon_certificate',false);
 save(fullfile(alt.output_dir,'measurement_report.mat'),'report');reports{i}=report;

end
if numel(reports)==1,reports=reports{1};end
end
