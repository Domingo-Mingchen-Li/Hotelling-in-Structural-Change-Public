function a=climate_parameter_audit(r,m)
s=r.trajectory;T=r.T;Cconv=zeros(1,T+1);
for t=0:T
 Cconv(t+1)=sum(m.rho.^t.*m.initial_components);
 for k=0:t-1,Cconv(t+1)=Cconv(t+1)+sum(m.a.*m.rho.^(t-1-k))*s.q(k+1);end
end
kernel=max(abs(Cconv-s.M))/max(1,max(abs(s.M)));
reconstruction=max(abs(s.C(:,m.L+1:end)-s.H(:,m.L+1:end)-s.future_components(:,m.L+1:end)),[],'all');
g=r.baseline.problem.growth;beta=r.baseline.p.beta;
tailerr=0;graderr=0;
for k=1:numel(m.etas)
 eta=m.etas(k);[v,d]=climate_environment_tail(s.C(:,end),s.H(:,end),s.q(end),g.g_e,beta,m,r.anchors.M_reference,eta);
 total=d.terminal_future_state+d.terminal_inherited_state;hist=d.terminal_inherited_state;num=0;
 for n=0:3000
  num=num+beta^n*((sum(total(1:4))/r.anchors.M_reference)^eta-(sum(hist(1:4))/r.anchors.M_reference)^eta);
  total=d.A*total;hist=d.A*hist;
 end
 tailerr=max(tailerr,abs(v-num)/max(1,abs(v)));
 h=1e-5*max(1,s.q(end));
 [vp,~]=climate_environment_tail(s.C(:,end),s.H(:,end),s.q(end)+h,g.g_e,beta,m,r.anchors.M_reference,eta);
 [vm,~]=climate_environment_tail(s.C(:,end),s.H(:,end),s.q(end)-h,g.g_e,beta,m,r.anchors.M_reference,eta);
 fd=(vp-vm)/(2*h);graderr=max(graderr,abs(fd-d.gradient_total_state(5))/max(1,abs(fd)));
end
gap=abs(r.baseline.validation.diagnostics.terminal.investment_growth_gap);
resource_mass=sum(r.baseline.path.R_input(m.L+1:T))+r.baseline.path.R_input(end)/(1-g.g_e);
mass_error=abs(resource_mass/r.inherited.R-1);
cal=max([r.scenarios.calibration_identity_error]);
support=max(cellfun(@(x)x.maximum_support_residual,r.welfare_audits));
a=struct('kernel_recursion_relative_error',kernel,'inherited_reconstruction_error',reconstruction, ...
 'environment_tail_relative_error',tailerr,'tail_gradient_error',graderr, ...
 'calibration_identity_error',cal,'supporting_price_residual',support,'baseline_terminal_gap',gap, ...
 'future_resource_mass_relative_error',mass_error);
a.passed=kernel<m.identity_tolerance&&reconstruction<m.identity_tolerance&& ...
 tailerr<m.identity_tolerance&&graderr<1e-7&&cal<m.identity_tolerance&&support<m.identity_tolerance&& ...
 gap<m.baseline_tail_gap_tolerance&&mass_error<1e-7;
assert(a.passed,'Climate:Audit','Environmental or welfare audit failed');

end
