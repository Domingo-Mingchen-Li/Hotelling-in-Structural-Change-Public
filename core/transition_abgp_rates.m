function g = transition_abgp_rates(p, mode)
a=[p.alpha_m,p.alpha_s,p.alpha_x,p.alpha_e];
b=[p.beta_m,p.beta_s,p.beta_x,p.beta_e]; gamma=1-a-b;
omega=[p.omega_m,p.omega_s]; epsv=p.epsilon;
LA=log([p.g_A_m,p.g_A_s,p.g_A_x,p.g_A_e]);
abar=omega*a(1:2)';bbar=omega*b(1:2)';
phi=(1/b(4))*(1-a(4)*b(3)/(1-a(3)))-epsv*(b(3)*abar/(1-a(3))+bbar);
g=struct('mode',mode,'legacy_phi',phi);
if strcmp(mode,'legacy')
    assert(abs(phi)>1e-12,'v2:ABGP','Legacy formula singular');
    pivot=(log(p.beta)+LA(4)/b(4)+epsv*LA(3)/(1-a(3))+epsv*(omega*LA(1:2)'))/phi;
    k=LA(3)/(1-a(3))+b(3)*pivot/(1-a(3));
    e=-LA(4)/b(4)-a(4)*LA(3)/((1-a(3))*b(4))+(1-a(4)*b(3)/(1-a(3)))*pivot/b(4);
    h=LA(4)/b(4)+(1+a(4)/b(4))*LA(3)/(1-a(3))+(b(3)*(1+a(4)/b(4))/(1-a(3))-1/b(4))*pivot;
    prices=-LA(1:2)+(1-a(1:2))*LA(3)/(1-a(3))+(b(3)*(1-a(1:2))/(1-a(3))-b(1:2))*pivot;
    P=omega*prices'; gross=-log(p.beta)+(1-epsv)*k+epsv*P;
    w=k; q=(LA(3)-gamma(3)*w)/b(3);
    g.linear_rcond=NaN;
elseif strcmp(mode,'paper')
    M=[gamma(3),b(3),0; -gamma(4),1,-b(4); ...
        -(1-epsv+epsv*(omega*gamma(1:2)')),-epsv*(omega*b(1:2)'),1];
    rhs=[LA(3);-LA(4);-log(p.beta)-epsv*(omega*LA(1:2)')];
    g.linear_rcond=rcond(M);
    assert(g.linear_rcond>1e-12,'v2:ABGP','Limiting equilibrium growth system singular');
    x=M\rhs;w=x(1);q=x(2);h=x(3);k=w;e=w-h;
    prices=-LA(1:2)+gamma(1:2)*w+b(1:2)*q;
    P=omega*prices';gross=h;g.linear_matrix=M;g.linear_rhs=rhs;
else
    error('v2:Mode','Unknown ABGP mode');
end
g.g_K=exp(k);g.g_E=exp(k);g.g_w=exp(w);g.g_q=exp(q);
g.g_e=exp(e);g.g_h=exp(h);g.g_pm=exp(prices(1));g.g_ps=exp(prices(2));g.g_P=exp(P);
g.r_ss=exp(gross)-(1-p.delta);
g.residual_names={'investment_unit_cost','processing_unit_cost','Euler','Hotelling','resource_revenue_growth'};
g.residuals=[gamma(3)*w+b(3)*q-LA(3), ...
    q+LA(4)-gamma(4)*w-b(4)*h, ...
    gross+log(p.beta)-(1-epsv)*k-epsv*P, h-gross, e+h-k];
g.nonhom_growth_paper=exp(-epsv*k+[p.zeta_m,p.zeta_s]*prices'+epsv*P);
g.nonhom_growth_legacy=exp(-epsv*k+[p.zeta_m,p.zeta_s]*prices'+2*epsv*P);
g.resource_tail_valid=isfinite(g.g_e)&&g.g_e>0&&g.g_e<1;
g.consistent=all(isfinite(g.residuals))&&max(abs(g.residuals))<1e-10&&g.resource_tail_valid;
end
