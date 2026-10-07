function report = calibrate_transition(c,resume_file)
if nargin<2,resume_file='';end
assert(strcmp(c.base.model,'paper'),'cal:Mode','Calibration requires paper mode');
assert(c.T>=50 && c.T==fix(c.T),'cal:Horizon','Invalid calibration horizon');
assert(isempty(c.check_T)||(isscalar(c.check_T)&&c.check_T>c.T&&c.check_T==fix(c.check_T)), ...
 'cal:Horizon','check_T must exceed T or be empty');
assert(isequal(size(c.targets),[3,1]) && all(isfinite(c.targets)) && ...
 all(c.scales>0) && numel(c.scales)==3 && all(c.fit_tolerances>0), ...
 'cal:Targets','Invalid targets, scales or tolerances');
assert(size(c.starts,2)==3 && all(isfinite(c.starts(:))) && ...
 all(c.lower<c.upper) && all(c.lower(2:3)>0) && ...
 all(all(c.starts>=c.lower & c.starts<=c.upper)), 'cal:Bounds','Invalid starts/bounds');
assert(c.fd_step>0&&c.minimum_step>0&&c.max_step>0&& ...
 c.initial_damping>0&&c.max_damping>=c.initial_damping&& ...
 c.continuation_initial_step>0&&c.continuation_initial_step<=1&& ...
 c.continuation_minimum_step>0&&c.continuation_minimum_step<=c.continuation_initial_step, ...
 'cal:Options','Invalid optimizer or continuation options');
base=c.base; base.horizons=c.T; base.policy_ids={}; base.policies=base.policies([]);
base.audit.derivatives=false; base.audit.restart=false; base=prepare_config(base);
raw=load(base.warm_start_file,'warm_starts');
seed=select_initial_path(raw.warm_starts,base.model,'baseline',c.T);
zero=struct('id','baseline','permanent',false,'start',0,'duration',11, ...
 'tau_c_m',0,'tau_c_s',0,'tau_int',0,'tau_e',0,'tau_x',0,'investment_policy_kind','none');
lower=encode(c.lower); upper=encode(c.upper);
cache={}; failures={}; evaluations=0;
if ~isfolder(c.output_root),mkdir(c.output_root);end
if isempty(resume_file)
 out=tempname(c.output_root);mkdir(out);
 save(fullfile(out,'calibration_configuration.mat'),'c');
else
 assert(isfile(resume_file),'cal:Resume','Checkpoint not found');
 out=fileparts(resume_file);
end
display_inner=isfield(c,'display_inner')&&logical(c.display_inner);
assert(isscalar(display_inner),'cal:Options','display_inner must be scalar');
logfid=fopen(fullfile(out,'calibration_live.txt'),'a');
assert(logfid>=0,'cal:Output','Cannot open incremental calibration progress log');
log_cleanup=onCleanup(@()fclose(logfid)); 
progress('Calibration progress log: %s\n',fullfile(out,'calibration_live.txt'));
report=struct('configuration',c,'output_dir',out,'starts',{{}}, ...
 'fit_passed',false,'identification_passed',false,'horizon_passed',false,'passed',false);
a=encode([seed.p.zeta_m,seed.Targets.K0,seed.Targets.R0]);
[ok,r,message]=inner(a,seed,c.T);
if ~ok,save(fullfile(out,'anchor_failure.mat'),'failures','message');end
assert(ok,'cal:Anchor','Warm start incompatible with fixed inputs: %s. Use validated current-growth baseline.',message);
cache{1}=struct('theta',a,'solution',r,'moments',moments(r.path));
best=[];
first_start=1;
if ~isempty(resume_file)
 checkpoint=load(resume_file,'report','best','failures');
 assert(isfield(checkpoint,'report')&&isfield(checkpoint,'best')&&~isempty(checkpoint.best), ...
  'cal:Resume','Checkpoint has no valid best equilibrium');
 assert(isequaln(checkpoint.report.configuration,c), ...
  'cal:Resume','Configuration changed. Restore the original config or start a fresh calibration.');
 report=checkpoint.report;best=checkpoint.best;failures=checkpoint.failures;
 report.resumed_from=resume_file;report.evaluation_count_since_resume=true;
 cache{end+1}=struct('theta',best.theta,'solution',best.solution,'moments',best.moments);
 first_start=numel(report.starts)+1;
 progress('Resumed %d completed starts from %s\n',numel(report.starts),resume_file);
end
for si=first_start:size(c.starts,1)
 progress('\nSTART %d/%d: constructing initial equilibrium\n',si,size(c.starts,1));
 theta=encode(c.starts(si,:)); [ok,m,r]=evaluate(theta);
 record=struct('start',c.starts(si,:),'status','initial_equilibrium_failed','history',[]);
 if ~ok
  report.starts{end+1}=record;save(fullfile(out,'progress.mat'),'report','best','failures');continue;
 end
 f=(m-c.targets)./c.scales; damping=c.initial_damping;
 history=[]; status='maximum_iterations';
 for it=0:c.max_iterations
  history(end+1,:)=[it,decode(theta),m.',sum(f.^2),damping]; 
  progress('start %d iter %d objective %.6g | zeta %.8g K0 %.8g R0 %.8g\n',si,it,sum(f.^2),decode(theta));
  if isempty(best)||sum(f.^2)<best.objective
   best=struct('theta',theta,'moments',m,'solution',r,'objective',sum(f.^2));
  end
  save(fullfile(out,'progress.mat'),'report','best','failures');
  if all(abs(m-c.targets)<=c.fit_tolerances),status='moment_tolerances_met';break;end
  if it==c.max_iterations,break;end
  [J,valid]=jacobian(theta,c.fd_step);
  if ~valid,status='jacobian_equilibrium_failed';break;end
  accepted=false;
  while damping<=c.max_damping
   d=-[J;sqrt(damping)*eye(3)]\[f;zeros(3,1)];
   d=d*min(1,c.max_step/max(norm(d,inf),eps));
   trial=min(upper,max(lower,theta+d));
   if norm(trial-theta,inf)<c.minimum_step,status='step_stalled';break;end
   [valid,mt,rt]=evaluate(trial);
   if valid
    ft=(mt-c.targets)./c.scales;
    if sum(ft.^2)<sum(f.^2)
     theta=trial;m=mt;r=rt;f=ft;damping=max(damping/3,1e-12);accepted=true;break;
    end
   end
   damping=damping*10;
  end
  if ~accepted
   if ~strcmp(status,'step_stalled'),status='damping_limit';end
   break;
  end
 end
 record.status=status;record.history=history;record.parameters=decode(theta);
 record.moments=m;record.fit_passed=all(abs(m-c.targets)<=c.fit_tolerances);
 report.starts{end+1}=record;
 save(fullfile(out,'progress.mat'),'report','best','failures');
end
assert(~isempty(best),'cal:NoSolution','No start produced a valid equilibrium. See configuration and warm start.');
report.parameters=decode(best.theta); report.moments=best.moments;
report.moment_errors=best.moments-c.targets;report.objective=best.objective;
report.fit_passed=all(abs(report.moment_errors)<=c.fit_tolerances);
progress('\nIDENTIFICATION AUDIT: finite-difference Jacobians\n');
[J1,v1]=jacobian(best.theta,c.fd_step);[J2,v2]=jacobian(best.theta,c.fd_step/2);
id=struct('passed',false,'coordinates','zeta_m, log(K0), log(R0)', ...
 'weighted_jacobian',J2,'step',c.fd_step/2);
if v1&&v2
 sv=svd(J2);id.singular_values=sv;id.condition_number=sv(1)/max(sv(end),realmin);
 id.relative_step_difference=norm(J1-J2,'fro')/max(norm(J2,'fro'),eps);
 id.rank=sum(sv>max(size(J2))*eps(max(sv)));
 id.at_bound=any(abs(best.theta-lower)<c.fd_step | abs(best.theta-upper)<c.fd_step);
 id.passed=id.rank==3 && id.condition_number<=c.maximum_condition_number && ...
  id.relative_step_difference<=c.jacobian_relative_tolerance && ~id.at_bound;
end
report.identification=id;report.identification_passed=id.passed;
matched=cellfun(@(s)isfield(s,'fit_passed')&&s.fit_passed,report.starts);
report.multistart=struct('checked',numel(report.starts)>1,'all_starts_fit',all(matched),'passed',false);
if all(matched)&&numel(report.starts)>1
 coordinates=cellfun(@(s)encode(s.parameters),report.starts,'UniformOutput',false);
 coordinates=cat(2,coordinates{:});
 report.multistart.maximum_spread=max(max(coordinates,[],2)-min(coordinates,[],2));
 report.multistart.passed=report.multistart.maximum_spread<=c.start_agreement_tolerance;
end
fitted=config_at(best.theta); fitted.audit.derivatives=true;fitted.audit.restart=true;
progress('\nFINAL EQUILIBRIUM AUDIT: derivatives and local restart\n');
audited=[]; 
log=evalc('audited=solve_transition(fitted,zero,c.T,best.solution,[]);');
report.equilibrium_audit=audited;report.audit_console=log;
report.finite_horizon_passed=audited.passed;
if audited.passed,best.solution=audited;end
report.horizon_checked=~isempty(c.check_T);
if report.horizon_checked
 progress('\nHORIZON AUDIT: T=%d\n',c.check_T);
 [valid,long,message]=inner(best.theta,best.solution,c.check_T);
 report.horizon=struct('T',c.check_T,'equilibrium_passed',valid,'error',message,'passed',false);
 if valid
  report.horizon.moments=moments(long.path);
  report.horizon.difference=report.horizon.moments-best.moments;
  report.horizon.passed=all(abs(report.horizon.difference)<=c.horizon_tolerances);
 end
 report.horizon_passed=report.horizon.passed;
end
report.passed=report.fit_passed&&report.identification_passed&&report.finite_horizon_passed ...
 && report.horizon_checked&&report.horizon_passed&&report.multistart.passed;
report.infinite_horizon_certificate=false;
report.evaluations=evaluations;report.failures=failures;
warm_starts={};
if best.solution.passed,warm_starts={compact(best.solution)};end
save(fullfile(out,'initial_paths.mat'),'warm_starts');
fitted.warm_start_file=fullfile(out,'initial_paths.mat');
save(fullfile(out,'fitted_configuration.mat'),'fitted');
save(fullfile(out,'calibration_report.mat'),'report');
fid=fopen(fullfile(out,'fitted_parameters.m'),'w');assert(fid>=0,'cal:Output','Cannot save parameter script');
fprintf(fid,'%% Apply to cfg=solver_config(); inspect calibration_report.mat before use.\n');
fprintf(fid,'cfg.p.zeta_m=%.17g; cfg.p.zeta_s=-cfg.p.zeta_m;\n',report.parameters(1));
fprintf(fid,'cfg.targets.K0=%.17g; cfg.targets.R0=%.17g;\n',report.parameters(2:3));
for sector={'m','s','x','e'}
 name=['g_A_' sector{1}];fprintf(fid,'cfg.p.%s=%.17g;\n',name,fitted.p.(name));
end
fprintf(fid,'cfg.warm_start_file=''%s'';\n',strrep(fitted.warm_start_file,'''',''''''));fclose(fid);
progress('Fit=%d identification=%d equilibrium audits=%d horizon=%d overall=%d\n', ...
 report.fit_passed,report.identification_passed,report.finite_horizon_passed,report.horizon_passed,report.passed);
progress('Results saved: %s\n',out);

 function progress(varargin)
  if startsWith(varargin{1},'Fit=')||startsWith(varargin{1},'Results saved:')
   fprintf(varargin{:});fprintf(logfid,varargin{:});
  end
  drawnow limitrate;
 end

 function cfg=config_at(t)
  p=decode(t);cfg=base;cfg.p.zeta_m=p(1);cfg.p.zeta_s=-p(1);
  cfg.targets.K0=p(2);cfg.targets.R0=p(3);cfg=prepare_config(cfg);
 end
 function [valid,s,msg]=inner(t,w,T)
  evaluations=evaluations+1;valid=false;s=[];msg='';
  progress('\nINNER %d BEGIN | T=%d zeta %.8g K0 %.8g R0 %.8g\n',evaluations,T,decode(t));
  clock=tic;
  try
   cfg=config_at(t);cfg.horizons=T;
   if display_inner
    cfg.newton_options.display=true;
    s=solve_transition(cfg,zero,T,w,[]);
   else
    captured=evalc('s=solve_transition(cfg,zero,T,w,[]);'); 
   end
   valid=s.passed;
   if ~valid,msg=s.error_message;end
  catch ME,msg=ME.message;end
  progress('INNER %d END | passed=%d elapsed=%.2fs\n',evaluations,valid,toc(clock));
  if ~valid,progress('Rejected inner attempt: %s\n',msg);end
  if ~valid,failures{end+1}=struct('parameters',decode(t),'T',T,'message',msg);end
 end
 function [valid,m,s]=evaluate(t)
  valid=false;m=nan(3,1);s=[];
  distances=cellfun(@(v)norm(v.theta-t),cache);[distance,k]=min(distances);
  if distance<1e-13,valid=true;m=cache{k}.moments;s=cache{k}.solution;return;end
  anchor=cache{k};[valid,s,~]=inner(t,anchor.solution,c.T);
  if ~valid
   progress('PARAMETER CONTINUATION: retry from nearest accepted path\n');
   lambda=0;increment=c.continuation_initial_step;w=anchor.solution;
   while lambda<1
    next=min(1,lambda+increment);tstage=anchor.theta+next*(t-anchor.theta);
    progress('Continuation lambda %.6g -> %.6g (step %.6g)\n',lambda,next,increment);
    [valid,st,~]=inner(tstage,w,c.T);
    if valid
     lambda=next;w=st;increment=min(2*increment,1-lambda);
    else
     increment=increment/2;
     progress('Continuation attempt rejected; next step %.6g\n',increment);
     if increment<c.continuation_minimum_step,break;end
    end
   end
   valid=lambda==1;if valid,s=w;end
  end
  if valid
   m=moments(s.path);cache{end+1}=struct('theta',t,'solution',s,'moments',m);
   if numel(cache)>60,cache(2)=[];end 
  end
 end
 function [J,valid]=jacobian(t,h)
  J=nan(3);valid=true;
  for j=1:3
   progress('JACOBIAN column %d/3: plus/minus parameter evaluations\n',j);
   tp=t;tm=t;tp(j)=min(upper(j),t(j)+h);tm(j)=max(lower(j),t(j)-h);
   [vp,mp]=evaluate(tp);[vm,mm]=evaluate(tm);
   if ~(vp&&vm)||tp(j)==tm(j),valid=false;return;end
   J(:,j)=((mp-mm)./c.scales)/(tp(j)-tm(j));
  end
 end
end
function t=encode(p),t=[p(1);log(p(2));log(p(3))];end
function p=decode(t),p=[t(1),exp(t(2)),exp(t(3))];end
function m=moments(p)
assert(numel(p.K)>=15,'cal:Moments','Need model dates 0:14');
lm=p.l_m([1,15]);ls=p.l_s([1,15]);den=lm+ls;
assert(all(isfinite([lm,ls]))&&all(lm>=0)&&all(ls>=0)&&all(den>0), ...
 'cal:Moments','Invalid labour inputs');
share=ls./den;Y=p.E(1)+p.I(1);
assert(isfinite(Y)&&Y>0,'cal:Moments','Invalid nominal final output');
m=[p.K(1)/Y;share(1);share(2)-share(1)];
assert(all(isfinite(m)),'cal:Moments','Nonfinite moments');
end
function s=compact(r)
s=struct('model_mode',r.model_mode,'case_name',r.case_name,'T',r.T,'p',r.p, ...
 'Targets',r.Targets,'policy',r.policy,'investment_policy_kind',r.investment_policy_kind, ...
 'boundary_kind',r.boundary_kind,'path',r.path,'validation',struct('passed',r.validation.passed));
end
