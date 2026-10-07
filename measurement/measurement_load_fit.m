function [cfg,entry,fit] = measurement_load_fit(m,id)
[cfg,~,~,signature,case_root]=measurement_case_config(m,id);
path=fullfile(case_root,'latest_calibration.mat');
assert(isfile(path),'measurement:MissingFit','Run run_measurement_robustness_calibration(''%s'') first',id);
raw=load(path,'entry');entry=raw.entry;
assert(entry.passed,'measurement:FitRejected','Latest calibration for %s is not accepted',id);
assert(isequaln(entry.signature,signature),'measurement:StaleFit', ...
 'Measurement inputs or calibration settings changed; recalibrate %s',id);
raw=load(entry.report_file,'report');fit=raw.report;
assert(fit.passed&&isequaln(fit.measurement_signature,signature), ...
 'measurement:FitRejected','Stored calibration report is not an accepted matching fit');
raw=load(entry.fitted_file,'fitted');fitted=raw.fitted;
for j=1:4
 sectors={'m','s','x','e'};
 for kind={'alpha','beta'}
  key=[kind{1} '_' sectors{j}];assert(fitted.p.(key)==cfg.p.(key), ...
   'measurement:StaleFit','Fitted production shares differ from current specification');
 end
end
cfg.p.zeta_m=fitted.p.zeta_m;cfg.p.zeta_s=-cfg.p.zeta_m;
fp=fitted.p;cp=cfg.p;
for name={'zeta_m','zeta_s','sect','A0'}
 if isfield(fp,name{1}),fp=rmfield(fp,name{1});end
 if isfield(cp,name{1}),cp=rmfield(cp,name{1});end
end
assert(isequaln(fp,cp)&&isequaln(fitted.A0,cfg.A0),'measurement:StaleFit', ...
 'Fitted fixed preferences/technology differ from current measurement configuration');
cfg.targets=fitted.targets;cfg.warm_start_file=entry.warm_start_file;
assert(isequal([cfg.p.zeta_m,cfg.targets.K0,cfg.targets.R0],fit.parameters), ...
 'measurement:FitRejected','Fitted parameters and report disagree');
assert(isfile(cfg.warm_start_file),'measurement:MissingFit','Fitted warm start missing');
cfg=prepare_config(cfg);
end
