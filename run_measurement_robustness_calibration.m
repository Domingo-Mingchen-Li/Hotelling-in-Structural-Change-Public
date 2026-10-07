function reports = run_measurement_robustness_calibration(case_id,resume_file,m)
setup_solver();
if nargin<3||isempty(m),m=measurement_robustness_config();end
if nargin<1||isempty(case_id),ids=m.case_ids;else,ids={char(case_id)};end
if nargin<2,resume_file='';else,resume_file=char(resume_file);end
assert(isempty(resume_file)||numel(ids)==1,'measurement:Resume','Resume one case at a time');
reports=cell(1,numel(ids));
for i=1:numel(ids)
 [~,c,spec,signature,case_root]=measurement_case_config(m,ids{i});
 if ~isfolder(c.output_root),mkdir(c.output_root);end
 entry=struct('passed',false,'signature',signature,'status','calibration_in_progress');
 save(fullfile(case_root,'latest_calibration.mat'),'entry');
 checkpoint=resume_file;
 if strcmp(checkpoint,'resume')
  files=dir(fullfile(c.output_root,'*','progress.mat'));
  assert(~isempty(files),'measurement:Resume','No checkpoint for %s',ids{i});
  [~,k]=max([files.datenum]);checkpoint=fullfile(files(k).folder,files(k).name);
 end
 if isempty(checkpoint),measurement_prepare_anchor(m,c,case_root);end
 assert(isfile(c.base.warm_start_file),'measurement:Anchor','Calibration anchor missing');

 report=calibrate_transition(c,checkpoint); 
 report.measurement=spec;report.measurement_signature=signature;
 save(fullfile(report.output_dir,'calibration_report.mat'),'report');
 raw=load(fullfile(report.output_dir,'fitted_configuration.mat'),'fitted');cfg=raw.fitted;
 fid=fopen(fullfile(report.output_dir,'measurement_fitted_parameters.m'),'w');
 assert(fid>=0,'measurement:Output','Cannot save complete measurement parameter script');
 fprintf(fid,'%% Complete production and fitted free parameters for measurement case %s.\n',spec.id);
 fprintf(fid,'%% Use only after calibration_report.passed is true. cfg=solver_config();\n');
 for sector={'m','s','x','e'}
  for kind={'alpha','beta'}
   key=[kind{1} '_' sector{1}];fprintf(fid,'cfg.p.%s=%.17g;\n',key,cfg.p.(key));
  end
 end
 fprintf(fid,'cfg.p.zeta_m=%.17g; cfg.p.zeta_s=-cfg.p.zeta_m;\n',cfg.p.zeta_m);
 fprintf(fid,'cfg.targets.K0=%.17g; cfg.targets.R0=%.17g;\n',cfg.targets.K0,cfg.targets.R0);
 fprintf(fid,'cfg.warm_start_file=''%s'';\n',strrep(cfg.warm_start_file,'''',''''''));fclose(fid);
 if isfile(fullfile(report.output_dir,'calibration_live.txt'))
  copyfile(fullfile(report.output_dir,'calibration_live.txt'), ...
   fullfile(report.output_dir,'measurement_console.txt'));
 end
 entry=struct('passed',report.passed,'signature',signature, ...
  'report_file',fullfile(report.output_dir,'calibration_report.mat'), ...
  'fitted_file',fullfile(report.output_dir,'fitted_configuration.mat'), ...
  'warm_start_file',fullfile(report.output_dir,'initial_paths.mat'));
 save(fullfile(case_root,'latest_calibration.mat'),'entry');reports{i}=report;
 assert(report.passed,'measurement:CalibrationFailed', ...
  'Fit/identification/multistart/equilibrium/horizon acceptance failed for %s; inspect %s',ids{i},entry.report_file);
end
if numel(reports)==1,reports=reports{1};end
end
