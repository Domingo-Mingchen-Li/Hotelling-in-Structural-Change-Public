function r=climate_pin_parameters(b,m,anchors)
P=b.path;p=b.p;T=b.T;g=b.problem.growth;
if nargin<3||isempty(anchors)
 e_ref=P.R_input(1);s=climate_trajectory(P.R_input,m,e_ref);
 Mref=s.M(m.reference_date+1);
 anchors=struct('e_reference',e_ref,'M_reference',Mref,'reference_date',m.reference_date);
else
 s=climate_trajectory(P.R_input,m,anchors.e_reference);Mref=anchors.M_reference;
end
assert(Mref>0&&isfinite(Mref),'Climate:Reference','Invalid burden reference');
rows=struct('eta',{},'loss_target',{},'chi',{},'discounted_damage',{},'utility_loss',{}, ...
 'environment_tail_share',{},'calibration_identity_error',{});
losses=zeros(size(m.loss_targets));wa=cell(size(m.loss_targets));
for j=1:numel(m.loss_targets)

 [losses(j),wa{j}]=climate_welfare_loss(b,m,m.loss_targets(j));

end
tail_details=cell(size(m.etas));
for k=1:numel(m.etas)
 eta=m.etas(k);ix=m.L+1:T;
 finite=sum(p.beta.^(0:T-m.L-1).*((s.M(ix)/Mref).^eta-(s.B(ix)/Mref).^eta));
 [tv,tail_details{k}]=climate_environment_tail(s.C(:,end),s.H(:,end),s.q(end),g.g_e,p.beta,m,Mref,eta);
 tail=p.beta^(T-m.L)*tv;H=finite+tail;
 assert(H>0&&isfinite(H),'Climate:Damage','Invalid discounted incremental damage');
 for j=1:numel(m.loss_targets)
  chi=losses(j)/H;err=abs(chi*H-losses(j))/max(losses(j),realmin);
  row=struct('eta',eta,'loss_target',m.loss_targets(j),'chi',chi,'discounted_damage',H, ...
   'utility_loss',losses(j),'environment_tail_share',tail/H,'calibration_identity_error',err);
  rows(end+1)=row; 

 end
end
idx=find([rows.eta]==m.benchmark_eta&[rows.loss_target]==m.benchmark_loss,1);
r=struct('T',T,'L',m.L,'anchors',anchors,'trajectory',s,'scenarios',rows,'benchmark',rows(idx), ...
 'welfare_audits',{wa},'environment_tail_details',{tail_details},'baseline',b, ...
 'economic_tail','geometric ABGP continuation from terminal E, prices and extraction');
r.inherited=struct('calendar_date',m.L,'K',P.K(m.L+1),'R',P.R(m.L+1), ...
 'A',[P.A_m(m.L+1);P.A_s(m.L+1);P.A_x(m.L+1);P.A_e(m.L+1)], ...
 'environment_components',s.inherited_components,'environment_burden',s.M(m.L+1));
end
