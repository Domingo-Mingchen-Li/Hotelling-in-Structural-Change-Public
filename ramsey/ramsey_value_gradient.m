function a=ramsey_value_gradient(x,q,ctx,need_gradient)
N=q.T;p=q.p;bb=p.beta;g=q.growth;X=reshape(x,5,N+1);
values=zeros(6,N+1);pols=cell(1,N+1);
for j=1:N+1
 pol=struct('tau_c_m',q.policy.tau_c_m_path(j),'tau_c_s',q.policy.tau_c_s_path(j), ...
  'tau_e',q.policy.tau_e_path(j),'tau_int',q.policy.tau_int_path(j),'tau_x',0,'investment_policy_kind','none');
 pols{j}=pol;values(:,j)=ramsey_local_objects(X(:,j),q.reference(:,j),q.technology(:,j),pol,p);
end
s=climate_trajectory(values(3,:),ctx.local_climate,ctx.anchors.e_reference);
eta=ctx.scenario.eta;chi=ctx.scenario.chi;Mref=ctx.anchors.M_reference;
weights=bb.^(0:N-1);gu=(g.g_E/g.g_P)^p.epsilon;gpi=g.g_pm^p.zeta_m*g.g_ps^p.zeta_s;
assert(bb*gu<1&&bb*gpi<1,'Ramsey:Tail','Private continuation diverges');
private=sum(weights.*(values(5,1:N)-values(6,1:N)))+bb^N*(values(5,end)/(1-bb*gu)-values(6,end)/(1-bb*gpi));
[tail,dt]=climate_environment_tail(s.C(:,end),zeros(4,1),s.q(end),g.g_e,bb,ctx.local_climate,Mref,eta);
damage=chi*(sum(weights.*(s.M(1:N)/Mref).^eta)+bb^N*tail);
a=struct('value',private-damage,'private_value',private,'discounted_damage',damage, ...
 'private_finite',sum(weights.*(values(5,1:N)-values(6,1:N))), ...
 'private_tail',bb^N*(values(5,end)/(1-bb*gu)-values(6,end)/(1-bb*gpi)), ...
 'damage_tail',chi*bb^N*tail,'environment',s,'local_values',values);
if ~need_gradient,return;end
pc=zeros(4,N+1);pc(:,end)=chi*bb^N*dt.gradient_total_state(1:4);
for t=N-1:-1:0
 pc(:,t+1)=chi*bb^t*eta*s.M(t+1)^(eta-1)/Mref^eta+ctx.local_climate.rho.*pc(:,t+2);
end
ew=[-ctx.local_climate.a'*pc(:,2:end)/ctx.anchors.e_reference, ...
 -chi*bb^N*dt.gradient_total_state(5)/ctx.anchors.e_reference];
aw=[weights,bb^N/(1-bb*gu)];bw=[-weights,-bb^N/(1-bb*gpi)];
WX=zeros(5,N+1);WP=zeros(N*numel(ctx.tools),1);FP=sparse(q.n,numel(WP));
h=ctx.options.derivative_step;
for j=1:N+1
 for k=1:5
  zp=X(:,j);zm=zp;zp(k)=zp(k)+h;zm(k)=zm(k)-h;
  dv=(ramsey_local_objects(zp,q.reference(:,j),q.technology(:,j),pols{j},p)- ...
      ramsey_local_objects(zm,q.reference(:,j),q.technology(:,j),pols{j},p))/(2*h);
  WX(k,j)=aw(j)*dv(5)+bw(j)*dv(6)+ew(j)*dv(3);
 end
 if j>N,continue;end 
 for k=1:numel(ctx.tools)
  name='tau_e';if strcmp(ctx.tools{k},'services'),name='tau_c_s';end
  pp=pols{j};pm=pp;pp.(name)=(1+pp.(name))*exp(h)-1;pm.(name)=(1+pm.(name))*exp(-h)-1;
  dv=(ramsey_local_objects(X(:,j),q.reference(:,j),q.technology(:,j),pp,p)- ...
      ramsey_local_objects(X(:,j),q.reference(:,j),q.technology(:,j),pm,p))/(2*h);
  col=(k-1)*N+j;WP(col)=aw(j)*dv(5)+bw(j)*dv(6)+ew(j)*dv(3);
  rows=5*(j-1)+(1:4);
  FP(rows,col)=[dv(1);-dv(2)/q.reference(1,j+1);dv(3)/q.reference(2,j);dv(4)];
  if j>1,FP(5*(j-2)+4,col)=-dv(4);end
 end
end
[J,~]=transition_sparse_jacobian(x,q,h);lambda=J'\WX(:);
a.adjoint_residual=norm(J'*lambda-WX(:),inf)/max(1,norm(WX(:),inf));
assert(a.adjoint_residual<1e-7,'Ramsey:Adjoint','Inaccurate adjoint linear solve');
a.gradient=(WP-FP'*lambda)/ctx.anchor;a.lambda=lambda;
a.environment_costate_discounted=pc;
end
