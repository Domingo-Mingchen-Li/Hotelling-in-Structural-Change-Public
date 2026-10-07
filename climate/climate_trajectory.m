function s=climate_trajectory(e,m,e_ref)
e=e(:).';assert(all(isfinite(e))&&all(e>=0)&&e_ref>0,'Climate:Flows','Invalid extraction');
T=numel(e)-1;assert(T>=m.L,'Climate:History','Missing history');
q=e/e_ref;C=zeros(4,T+1);C(:,1)=m.initial_components;
for j=1:T,C(:,j+1)=m.rho.*C(:,j)+m.a*q(j);end
CL=C(:,m.L+1);H=nan(4,T+1);
for t=m.L:T,H(:,t+1)=m.rho.^(t-m.L).*CL;end
s=struct('dates',0:T,'q',q,'C',C,'M',sum(C,1),'H',H,'B',sum(H,1), ...
 'inherited_components',CL,'e_reference',e_ref);
s.future_components=C-H;
end
