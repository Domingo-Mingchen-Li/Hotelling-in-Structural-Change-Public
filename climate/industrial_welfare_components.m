function a=industrial_welfare_components(policy,baseline,m,anchors,scenarios,w)
assert(policy.T==baseline.T,'Welfare:Horizon','Mismatched local horizons');
p=baseline.p;g=baseline.problem.growth;N=baseline.T;bb=p.beta;
[uP,tailP]=private_value(policy);[uB,tailB]=private_value(baseline);
weights=bb.^(0:N-1);du=uP-uB;
private_finite=sum(weights.*du(1:N));private_tail=bb^N*(tailP-tailB);
private=private_finite+private_tail;
sP=climate_trajectory(policy.path.R_input,m,anchors.e_reference);
sB=climate_trajectory(baseline.path.R_input,m,anchors.e_reference);
Mref=anchors.M_reference;
anchor=scenarios(find([scenarios.eta]==m.benchmark_eta&[scenarios.loss_target]==m.benchmark_loss,1)).utility_loss;
assert(isfinite(anchor)&&anchor>0,'Welfare:Anchor','Missing positive utility-loss anchor');
rows=[];envpaths=zeros(numel(scenarios),N+1);tailchecks=zeros(1,numel(scenarios));
for k=1:numel(scenarios)
 sc=scenarios(k);eta=sc.eta;chi=sc.chi;
 DD=(sP.M/Mref).^eta-(sB.M/Mref).^eta;
 [tp,dp]=climate_environment_tail(sP.C(:,end),sP.H(:,end),sP.q(end),g.g_e,bb,m,Mref,eta);
 [tb,db]=climate_environment_tail(sB.C(:,end),sB.H(:,end),sB.q(end),g.g_e,bb,m,Mref,eta);
 if eta==2
  yp=[sP.C(:,end);sP.q(end)];yb=[sB.C(:,end);sB.q(end)];
  direct=(yp'*dp.P*yp-yb'*db.P*yb)/Mref^2;
 else
  direct=dp.linear_value'*([sP.C(:,end);sP.q(end)]-[sB.C(:,end);sB.q(end)])/Mref;
 end
 tailchecks(k)=abs((tp-tb)-direct)/max(1,abs(direct));
 env_finite=-chi*sum(weights.*DD(1:N));env_tail=-chi*bb^N*(tp-tb);
 env=env_finite+env_tail;total=private+env;
 row=struct('policy',policy.case_name,'calendar_terminal',N+m.announcement_calendar_date, ...
  'announcement_date',m.announcement_calendar_date,'eta',eta,'loss_target',sc.loss_target,'chi',chi, ...
  'consumption_utility_change',private,'environment_welfare_change',env,'net_welfare_change',total, ...
  'consumption_finite',private_finite,'consumption_tail',private_tail, ...
  'environment_finite',env_finite,'environment_tail',env_tail, ...
  'consumption_in_anchor_units',private/anchor,'environment_in_anchor_units',env/anchor, ...
  'net_in_anchor_units',total/anchor,'consumption_sign',sgn(private,anchor,w), ...
  'environment_sign',sgn(env,anchor,w),'net_sign',sgn(total,anchor,w), ...
  'finite_cumulative_extraction_percent',100*(sum(policy.path.R_input(1:N))/sum(baseline.path.R_input(1:N))-1), ...
  'remaining_resource_difference',policy.path.R(end)-baseline.path.R(end), ...
  'tail_identity_error',tailchecks(k));
 if isempty(rows),rows=row;else,rows(end+1)=row;end 
 envpaths(k,:)=-chi*DD;

end
assert(max(tailchecks)<1e-9,'Welfare:Tail','Environmental tail identity failed');
assert(all(isfinite([private;envpaths(:)])),'Welfare:Finite','Nonfinite welfare');
a=struct('rows',rows,'private_flow_difference',du,'environment_flow_difference',envpaths, ...
 'policy_environment',sP,'baseline_environment',sB,'passed',true, ...
 'unit_note','Utility differences and multiples of the FIXED benchmark basket-loss anchor; NOT CE percentages', ...
 'tail_note','Exact discounted continuation conditional on constant terminal wedges and geometric ABGP prices, E and extraction');
end
function [u,tail]=private_value(r)
p=r.p;g=r.problem.growth;P=r.path;
pm=P.p_m.*(1+r.policy.tau_c_m_path);ps=P.p_s.*(1+r.policy.tau_c_s_path);
lP=p.omega_m*log(pm)+p.omega_s*log(ps);lPi=p.zeta_m*log(pm)+p.zeta_s*log(ps);
a=exp(p.epsilon*(log(P.E)-lP))/p.epsilon;b=p.xi*exp(lPi);u=a-b;
gu=(g.g_E/g.g_P)^p.epsilon;gpi=g.g_pm^p.zeta_m*g.g_ps^p.zeta_s;
assert(p.beta*gu<1&&p.beta*gpi<1,'Welfare:Tail','Private welfare diverges');
tail=a(end)/(1-p.beta*gu)-b(end)/(1-p.beta*gpi);
err=max(abs((pm.*P.c_m+ps.*P.c_s)./P.E-1));
assert(err<1e-8,'Welfare:Prices','Purchaser-price expenditure identity failed');
end
function v=sgn(x,anchor,w)
v=sign(x);if abs(x)<=w.sign_anchor_tolerance*anchor,v=0;end
end
