function [N,audit]=climate_welfare_loss(b,m,lambda)
p=b.p;g=b.problem.growth;P=b.path;T=b.T;
if strcmp(m.loss_definition,'fixed_price_expenditure')
 fac=(1-(1-lambda)^p.epsilon)/p.epsilon;
 u=(P.E./P.P_star).^p.epsilon;
 ratio=p.beta*(g.g_E/g.g_P)^p.epsilon;
 N=fac*(sum(p.beta.^(0:T-m.L-1).*u(m.L+1:T))+ ...
  p.beta^(T-m.L)*u(T+1)/(1-ratio));
 audit=struct('definition',m.loss_definition,'tail_years',Inf,'tail_change_relative',0, ...
  'maximum_support_residual',0,'tail_exact',true,'tail_numerically_converged',true);
 return;
end
N=0;maxres=0;
for t=m.L:T-1
 [loss,a]=climate_basket_loss(P.c_m(t+1),P.c_s(t+1),P.p_m(t+1),P.p_s(t+1),p,lambda);
 N=N+p.beta^(t-m.L)*loss;maxres=max(maxres,a.residual);
end
baseN=N;tail=0;last=0;smallblocks=0;change=Inf;converged=false;
for k=0:m.maximum_utility_tail_years-1
 pm=P.p_m(end)*g.g_pm^k;ps=P.p_s(end)*g.g_ps^k;E=P.E(end)*g.g_E^k;
 PP=exp(p.omega_m*log(pm)+p.omega_s*log(ps));Pi=exp(p.zeta_m*log(pm)+p.zeta_s*log(ps));
 term=p.xi*E^(1-p.epsilon)*PP^p.epsilon*Pi;
 cm=(p.omega_m*E+p.zeta_m*term)/pm;cs=(p.omega_s*E+p.zeta_s*term)/ps;
 assert(cm>0&&cs>0&&all(isfinite([cm cs pm ps E])),'Climate:Tail','Invalid utility-tail basket');
 [loss,a]=climate_basket_loss(cm,cs,pm,ps,p,lambda);
 maxres=max(maxres,a.residual);tail=tail+p.beta^(T-m.L+k)*loss;
 if mod(k+1,m.utility_tail_block_years)==0
  change=abs(tail-last)/max(abs(baseN+tail),realmin);last=tail;
  if change<m.utility_tail_relative_tolerance,smallblocks=smallblocks+1;else,smallblocks=0;end
  if smallblocks>=2,converged=true;break;end
 end
end
assert(converged,'Climate:UtilityTail','Utility-loss tail not converged; raise maximum_utility_tail_years');
N=baseN+tail;
audit=struct('definition',m.loss_definition,'tail_years',k+1,'tail_change_relative',change, ...
 'maximum_support_residual',maxres,'tail_exact',false,'tail_numerically_converged',converged, ...
 'tail_method','two successive block increments below tolerance; numerical convergence, not rigorous bound');
end
