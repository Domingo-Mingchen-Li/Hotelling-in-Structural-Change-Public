function report = run_solver(cfg)
cfg=prepare_config(cfg);assert(isfile(cfg.warm_start_file),'v2:WarmStart','Warm start file not found');
raw=load(cfg.warm_start_file,'warm_starts');pool=raw.warm_starts(:).';
[out,started]=create_run_directory(cfg);
transcript=evalc('report=run_cases(cfg,pool);');report.output_dir=out;
report.configuration=cfg;report.run_started_at=started;
if isfield(cfg,'run_kind'),report.run_kind=cfg.run_kind;else,report.run_kind='experiments';end
fid=fopen(fullfile(out,'console.txt'),'w','n','UTF-8');assert(fid>=0,'v2:Output','Cannot save console');
fwrite(fid,transcript,'char');fclose(fid);
save(fullfile(out,'configuration.mat'),'cfg');save(fullfile(out,'report.mat'),'report');
baselines=report.baseline_cases;solutions={};warm_starts={};
for j=1:numel(baselines),if baselines{j}.passed,warm_starts{end+1}=compact_path(baselines{j});end,end 
for k=1:numel(report.policy_cases)
    for j=1:numel(report.policy_cases{k}.results)
        s=report.policy_cases{k}.results{j};solutions{end+1}=s; 
        if s.passed,warm_starts{end+1}=compact_path(s);end 
    end
end
save(fullfile(out,'solutions.mat'),'baselines','solutions');
save(fullfile(out,'initial_paths.mat'),'warm_starts');
if report.passed
    try
        export_results(report,out);
    catch ME
        report.export_error=ME.message;save(fullfile(out,'report.mat'),'report');

    end
end

if ~report.passed,error('v2:RunFailed','One or more cases failed; see saved console/report');end
end

function report=run_cases(cfg,pool)
report=struct('passed',false,'model_mode',cfg.model,'horizons',cfg.horizons, ...
    'baseline_cases',{{}},'policy_cases',{{}},'empirical_fit_validated',false, ...
    'infinite_horizon_accuracy_validated',false);
zero=struct('id','baseline','permanent',false,'start',0,'duration',11, ...
    'tau_c_m',0,'tau_c_s',0,'tau_int',0,'tau_e',0,'tau_x',0,'investment_policy_kind','none');
for j=1:numel(cfg.horizons)
    T=cfg.horizons(j);report.baseline_cases{j}=select_and_solve(cfg,zero,T,pool,[]);
end
if ~all(cellfun(@(s)s.passed,report.baseline_cases)),return;end
for k=1:numel(cfg.policies)
    spec=cfg.policies(k);record=struct('id',spec.id,'specification',spec,'results',{{}});
    for j=1:numel(cfg.horizons)
        T=cfg.horizons(j);record.results{j}=select_and_solve(cfg,spec,T,pool,report.baseline_cases{j});
    end
    record.passed=all(cellfun(@(s)s.passed,record.results));report.policy_cases{k}=record;
end
report.passed=all(cellfun(@(s)s.passed,report.policy_cases));
if report.passed&&cfg.audit.horizons
    try,report=audit_horizons(report);
    catch ME
        report.passed=false;report.audit_error=getReport(ME,'extended','hyperlinks','off');

    end
end
if report.passed

end
end

function s=compact_path(r)
s=struct('model_mode',r.model_mode,'case_name',r.case_name,'T',r.T,'p',r.p, ...
    'Targets',r.Targets,'policy',r.policy,'investment_policy_kind',r.investment_policy_kind, ...
    'boundary_kind',r.boundary_kind,'path',r.path,'validation',struct('passed',r.validation.passed));
end

function r=select_and_solve(cfg,spec,T,pool,baseline)
try
    warm=select_initial_path(pool,cfg.model,spec.id,T);
    r=solve_transition(cfg,spec,T,warm,baseline);
catch ME
    r=struct('passed',false,'case_name',spec.id,'model_mode',cfg.model,'T',T, ...
        'error_identifier',ME.identifier,'error_message',ME.message);

end
end
