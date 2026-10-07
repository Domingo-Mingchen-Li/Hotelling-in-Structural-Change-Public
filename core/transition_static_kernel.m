function s = transition_static_kernel(z, tech, pol, p, mode)
assert(any(strcmp(mode,{'legacy','paper'})), 'v2:Mode','Unknown model mode');
for name={'K','R','r','h','E'}
    x=z.(name{1});
    assert(isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0, ...
        'v2:State','Invalid state %s',name{1});
end
for name={'tau_c_m','tau_c_s','tau_int','tau_e','tau_x'}
    x=pol.(name{1});
    assert(isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>-1, ...
        'v2:Policy','Invalid wedge %s',name{1});
end
purchase=isfield(pol,'investment_policy_kind')&&strcmp(pol.investment_policy_kind,'purchase_subsidy');
assert(~purchase||strcmp(mode,'paper'),'v2:Scope','Purchase subsidy requires paper mode');
assert(strcmp(mode,'legacy')||pol.tau_x==0||purchase,'v2:Scope', ...
    'Nonzero paper tau_x requires an explicit investment purchase specification');
qI=1; if purchase,qI=1+pol.tau_x;end
sectors={'m','s','x','e'}; a=zeros(1,4); b=a; gamma=a; A=a; d=a;
for j=1:4
    i=sectors{j}; a(j)=p.(['alpha_' i]); b(j)=p.(['beta_' i]);
    gamma(j)=1-a(j)-b(j); A(j)=tech.(['A_' i]);
    assert(a(j)>0&&b(j)>0&&gamma(j)>0&&A(j)>0&&isfinite(A(j)), ...
        'v2:Technology','Invalid production coefficients or technology');
    d(j)=a(j)*log(a(j))+b(j)*log(b(j))+gamma(j)*log(gamma(j));
end
lr=log(z.r); lh=log(z.h)+log1p(pol.tau_e); wedge=log1p(pol.tau_int);
den=gamma(3)+b(3)*gamma(4);
lw=(log(A(3))+d(3)-a(3)*lr-b(3)*(a(4)*lr+b(4)*lh-log(A(4))-d(4)+wedge))/den;
lq=a(4)*lr+gamma(4)*lw+b(4)*lh-log(A(4))-d(4);
lprices=a(1:3)*lr+gamma(1:3)*lw+b(1:3)*(lq+wedge)-log(A(1:3))-d(1:3);
s=struct(); s.mode=mode; s.w=exp(lw); s.p_int=exp(lq);
s.p_int_user=exp(lq+wedge); s.p_m=exp(lprices(1)); s.p_s=exp(lprices(2));
lg=[lprices(1)+log1p(pol.tau_c_m),lprices(2)+log1p(pol.tau_c_s)];
omega=[p.omega_m,p.omega_s]; zeta=[p.zeta_m,p.zeta_s];
logP=sum(omega.*lg); logPi=sum(zeta.*lg);
s.P_star=exp(logP); s.Pi_star_paper=exp(logPi);
power=1; if strcmp(mode,'legacy'), power=2; end
term=p.xi*exp((1-p.epsilon)*log(z.E)+logPi+power*p.epsilon*logP);
raw_c=(omega*z.E+zeta*term)./exp(lg);
s.raw_consumption=raw_c; s.clamp_active=false;
if strcmp(mode,'legacy')
    s.clamp_active=any(raw_c<1e-9); c=max(raw_c,1e-9);
else
    c=raw_c;
end
s.c_m=c(1); s.c_s=c(2); s.shares=exp(lg).*c/z.E;
Vcons=[s.p_m*s.c_m,s.p_s*s.c_s]; pass=1/(1+pol.tau_int);
s.I=(s.w-sum((gamma(1:2)+gamma(4)*b(1:2)*pass).*Vcons))/ ...
    (gamma(3)+gamma(4)*b(3)*pass);
Vfinal=[Vcons,s.I]; Ve=sum(b(1:3).*Vfinal)*pass;
V=[Vfinal,Ve]; s.k=a.*V/z.r; s.l=gamma.*V/s.w;
s.m_inputs=b(1:3).*Vfinal/s.p_int_user; s.m_total=Ve/s.p_int;
s.e=b(4)*Ve/(z.h*(1+pol.tau_e)); s.K_d=sum(s.k);
s.y=[c,s.I,s.m_total]; s.V=V;
s.Y_final=sum(Vfinal);
s.transfer=pol.tau_c_m*Vcons(1)+pol.tau_c_s*Vcons(2)+ ...
    pol.tau_int*s.p_int*sum(s.m_inputs)+pol.tau_e*z.h*s.e;
s.investment_user_price=qI;
s.investment_transfer=0;
if purchase,s.investment_transfer=pol.tau_x*s.I;end
s.transfer=s.transfer+s.investment_transfer;
s.budget_raw=s.w+z.r*z.K+z.h*s.e+s.transfer-qI*s.I-z.E;
s.feasible=all(isfinite([s.w,s.p_int,s.p_m,s.p_s,s.P_star,c,s.I,s.e,s.k,s.l,s.m_inputs])) ...
    && all(c>0)&&s.I>0&&s.e>0&&all(s.k>0)&&all(s.l>0)&&all(s.m_inputs>0) ...
    &&s.e<=z.R&&~s.clamp_active;
s.residual_names={'labor','intermediate','expenditure','budget','numeraire','production'};
prod_res=NaN;
if all(s.y>0)&&all(s.k>0)&&all(s.l>0)&&all(s.m_inputs>0)&&s.e>0
    inputs=[s.m_inputs,s.e];
    prod_res=max(abs(log(s.y)-log(A)-a.*log(s.k)-gamma.*log(s.l)-b.*log(inputs)));
end
s.residuals=[sum(s.l)-1, (sum(s.m_inputs)-s.m_total)/max(abs(s.m_total),eps), ...
    (sum(exp(lg).*c)-z.E)/z.E, s.budget_raw/max(abs(z.E)+abs(s.I),eps), ...
    lprices(3),prod_res];
s.capital_residual=(s.K_d-z.K)/z.K;
end
