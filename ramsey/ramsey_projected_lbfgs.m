function o=ramsey_projected_lbfgs(u,ctx,warm,checkpoint)
w=ctx.options;lo=ctx.lower;hi=ctx.upper;u=min(hi,max(lo,u(:)));
e=ramsey_evaluate(u,ctx,warm,true);assert(e.success,'Ramsey:Seed','Initial policy failed');
S={};Y={};history=[];status='iteration_limit';passed=false;
for it=0:w.max_iterations
 g=e.gradient;pg=norm(u-min(hi,max(lo,u-g)),inf);
 row=struct('iteration',it,'gain_in_anchor_units',e.gain,'projected_gradient',pg, ...
  'equilibrium_residual',e.maximum_residual,'max_resource_tax',max(e.problem.policy.tau_e_path), ...
  'max_abs_service_wedge',max(abs(e.problem.policy.tau_c_s_path)));
 if isempty(history),history=row;else,history(end+1)=row;end 

 snapshot=struct('evaluation',e,'history',history,'status','running','passed',false);
 save(checkpoint,'snapshot');
 if pg<=w.projected_gradient_tolerance,passed=true;status='projected_stationarity';break;end
 if it==w.max_iterations,break;end
 r=g;al=zeros(1,numel(S));
 for k=numel(S):-1:1,al(k)=(S{k}'*r)/(S{k}'*Y{k});r=r-al(k)*Y{k};end
 scale=1;if ~isempty(S),scale=min(1e3,max(1e-3,(S{end}'*Y{end})/(Y{end}'*Y{end})));end
 r=scale*r;
 for k=1:numel(S),be=(Y{k}'*r)/(S{k}'*Y{k});r=r+S{k}*(al(k)-be);end
 d=-r;d(u<=lo+1e-12&d<0)=0;d(u>=hi-1e-12&d>0)=0;
 if g'*d>=-1e-12*norm(g)*max(1,norm(d)),d=min(hi,max(lo,u-g))-u;end
 accepted=false;
 for attempt=1:2
  if attempt==2

   S={};Y={};d=min(hi,max(lo,u-g))-u;
  end
  step=min(1,w.maximum_control_step/max(norm(d,inf),eps));
  for bt=0:w.max_backtracks
   un=min(hi,max(lo,u+step*d));move=un-u;
   if norm(move,inf)<1e-12,break;end
   en=ramsey_evaluate(un,ctx,e.x,true);
   if en.success&&en.f<=e.f+w.armijo*(g'*move),accepted=true;break;end
   step=step/2;
  end
  if accepted,break;end
 end
 if ~accepted,status='line_search_failed';break;end
 ss=un-u;yy=en.gradient-g;
 if ss'*yy>1e-12*norm(ss)*norm(yy)
  S{end+1}=ss;Y{end+1}=yy; 
  if numel(S)>w.memory,S=S(2:end);Y=Y(2:end);end
 end
 u=un;e=en;
end
o=struct('passed',passed,'status',status,'evaluation',e,'history',history, ...
 'projected_gradient',pg,'iterations',it,'global_optimum_certified',false);
snapshot=o;save(checkpoint,'snapshot');
end
