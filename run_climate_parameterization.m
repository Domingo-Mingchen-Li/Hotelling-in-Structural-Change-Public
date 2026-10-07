function report=run_climate_parameterization(m)
setup_solver();if nargin<1||isempty(m),m=climate_config();end
m=climate_prepare_config(m);cfg=m.base;
cfg.warm_start_file=m.seed_file;cfg.output_root=m.output_root;
cfg.run_kind='experiments';cfg.output_tag='climate_parameters';cfg.make_plots=false;
cfg.newton_options.display=m.display_inner;cfg.audit.derivatives=m.audit_derivatives;
cfg.audit.restart=m.audit_restart;[out,started]=create_run_directory(cfg);
report=struct('passed',false,'configuration',m,'output_dir',out,'started_at',started, ...
 'policy_optimized',false,'physical_emissions_calibrated',false);
entry=struct('passed',false,'report_file',fullfile(out,'climate_report.mat'), ...
 'parameter_file',fullfile(out,'climate_parameters.mat'));
save(fullfile(m.output_root,'latest_climate.mat'),'entry');
previous=get(0,'Diary');previous_file=get(0,'DiaryFile');
diary(fullfile(out,'console.txt'));diary on;
cleanup=onCleanup(@()restore_diary(previous,previous_file)); 

try
 assert(isfile(m.seed_file),'Climate:Seed','Missing seed file: %s',m.seed_file);
 raw=load(m.seed_file,'warm_starts');pool=raw.warm_starts(:).';
 spec=struct('id','baseline','start',0,'duration',10,'permanent',false, ...
  'tau_c_m',0,'tau_c_s',0,'tau_e',0,'tau_int',0,'tau_x',0,'investment_policy_kind','none');

 warm=select_initial_path(pool,'paper','baseline',m.T);
 b=solve_transition(cfg,spec,m.T,warm,[]);assert(b.passed,'Climate:Baseline','Main baseline failed');
 report.main_baseline=b;save(fullfile(out,'climate_report.mat'),'report');

 c=solve_transition(cfg,spec,m.check_T,b,[]);assert(c.passed,'Climate:Baseline','Check baseline failed');
 report.check_baseline=c;save(fullfile(out,'climate_report.mat'),'report');

 report.main=climate_pin_parameters(b,m,[]);
 report.main_audit=climate_parameter_audit(report.main,m);save(fullfile(out,'climate_report.mat'),'report');

 report.check=climate_pin_parameters(c,m,report.main.anchors);
 report.check_audit=climate_parameter_audit(report.check,m);
 report.horizon_audit=climate_horizon_audit(report.main,report.check,m);
 report.passed=report.main_audit.passed&&report.check_audit.passed&&report.horizon_audit.passed;
 report.completed_at=char(datetime('now'));

 report.exports=climate_export_parameters(report,out);
 save(fullfile(out,'climate_report.mat'),'report');entry.passed=report.passed;
 save(fullfile(m.output_root,'latest_climate.mat'),'entry');

 assert(report.passed,'Climate:Acceptance','Horizon check failed; saved diagnostics are not accepted parameters');
catch ME
 report.passed=false;report.error_identifier=ME.identifier;report.error_message=ME.message;
 save(fullfile(out,'climate_report.mat'),'report');entry.passed=false;
 save(fullfile(m.output_root,'latest_climate.mat'),'entry');
 rethrow(ME);
end
end
function restore_diary(state,file)
diary off;if strcmpi(state,'on'),diary(file);diary on;end
end
