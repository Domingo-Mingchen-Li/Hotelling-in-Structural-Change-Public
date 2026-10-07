function a=ramsey_gradient_audit(e,ctx)
n=numel(e.u);rows=[];h=ctx.options.gradient_audit_step;
for k=1:2
 direction=sin((1:n)'*(0.71+0.13*k));direction=direction/max(abs(direction));
 ep=ramsey_evaluate(e.u+h*direction,ctx,e.x,false);
 em=ramsey_evaluate(e.u-h*direction,ctx,e.x,false);
 assert(ep.success&&em.success,'Ramsey:GradientAudit','Directional equilibrium probe failed');
 fd=(ep.f-em.f)/(2*h);ad=e.gradient'*direction;
 err=abs(fd-ad)/max([1 abs(fd) abs(ad)]);
 row=struct('direction',k,'finite_difference',fd,'adjoint',ad,'error',err, ...
  'passed',err<ctx.options.gradient_audit_tolerance);
 if isempty(rows),rows=row;else,rows(end+1)=row;end 
end
a=struct('passed',all([rows.passed]),'directions',rows,'maximum_error',max([rows.error]));

assert(a.passed,'Ramsey:GradientAudit','Economic/environmental adjoint failed');
end
