function e=ramsey_evaluate(u,ctx,warm,need_gradient)
if nargin<4,need_gradient=true;end
q=ctx.problem;N=q.T;assert(numel(u)==N*numel(ctx.tools),'Ramsey:Controls','Incorrect control dimension');
for k=1:numel(ctx.tools)
 name='tau_e_path';if strcmp(ctx.tools{k},'services'),name='tau_c_s_path';end
 q.policy.(name)=[expm1(u((k-1)*N+(1:N))).',0];
end
e=struct('success',false,'u',u(:),'x',warm,'problem',q);
try
 opt=ctx.options.inner_options;opt.display=ctx.options.display_inner;
 sol=transition_sparse_newton(q,warm,opt);
 assert(sol.success,'Ramsey:Equilibrium','Inner equilibrium: %s',sol.status);
 [F,d]=transition_global_residual(sol.x,q);
 assert(d.all_feasible&&norm(F,inf)<=1e-9,'Ramsey:Equilibrium','Feasibility or residual check failed');
 a=ramsey_value_gradient(sol.x,q,ctx,need_gradient);
 e.success=true;e.x=sol.x;e.solver=sol;e.welfare=a;
 e.gain=(a.value-ctx.baseline_value)/ctx.anchor;e.f=-e.gain;
 if need_gradient,e.gradient=-a.gradient;end
 e.maximum_residual=norm(F,inf);e.terminal_gap=d.terminal.investment_growth_gap;

catch ME
 e.error_identifier=ME.identifier;e.error_message=ME.message;

end
end
