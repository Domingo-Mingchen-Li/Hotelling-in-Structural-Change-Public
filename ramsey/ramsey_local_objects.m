function [v,s]=ramsey_local_objects(logz,reference,technology,pol,p)
z=reference(:).*exp(logz(:));
st=struct('K',z(1),'R',z(2),'r',z(3),'h',z(4),'E',z(5));
tech=struct('A_m',technology(1),'A_s',technology(2),'A_x',technology(3),'A_e',technology(4));
s=transition_static_kernel(st,tech,pol,p,'paper');
lpm=log(s.p_m)+log1p(pol.tau_c_m);lps=log(s.p_s)+log1p(pol.tau_c_s);
lP=p.omega_m*lpm+p.omega_s*lps;lPi=p.zeta_m*lpm+p.zeta_s*lps;
a=exp(p.epsilon*(log(z(5))-lP))/p.epsilon;b=p.xi*exp(lPi);
mu=(p.epsilon-1)*log(z(5))-p.epsilon*lP;
v=[s.capital_residual;s.I;s.e;mu;a;b];
end
