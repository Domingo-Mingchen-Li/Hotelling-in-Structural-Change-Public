function a = audit_policy_equations(path,q)
p=q.p;pol=q.policy;t=1:q.T;n=t+1;
mu=path.E.^(p.epsilon-1)./path.P_star.^p.epsilon;
if strcmp(q.mode,'legacy'),mu=mu./(1+pol.tau_x_path);end
qI=ones(size(q.time));
if strcmp(q.investment_policy_kind,'purchase_subsidy'),qI=1+pol.tau_x_path;end
payoff=path.r(n)+(1-p.delta)*qI(n);
left=qI(t).*mu(t);right=p.beta*mu(n).*payoff;
euler=(left-right)./max(abs(left)+abs(right),realmin);
left=mu(t).*path.h(t);right=p.beta*mu(n).*path.h(n);
resource_foc=(left-right)./max(abs(left)+abs(right),realmin);
V=[path.p_m.*path.c_m;path.p_s.*path.c_s;path.I];
Ve=([p.beta_m,p.beta_s,p.beta_x]*V)./(1+pol.tau_int_path);
Tcons=pol.tau_c_m_path.*V(1,:)+pol.tau_c_s_path.*V(2,:);
Tint=pol.tau_int_path.*Ve;Te=pol.tau_e_path.*path.h.*path.R_input;
Tx=(qI-1).*path.I;transfer=Tcons+Tint+Te+Tx;
income=path.w+path.r.*path.K+path.h.*path.R_input;
budget=(income+transfer-path.E-qI.*path.I)./(path.E+abs(qI.*path.I));
transfer_error=(path.transfer-transfer)./(path.E+abs(path.I));
net_resources=(income+Tcons+Tint+Te-path.E-path.I)./(path.E+abs(path.I));
expected=payoff./qI(t);rent=(path.h(n)./path.h(t)-expected)./expected;
a=struct('passed',false,'maximum_euler_foc',max(abs(euler)), ...
    'maximum_resource_foc',max(abs(resource_foc)),'maximum_rent_return_error',max(abs(rent)), ...
    'maximum_household_budget',max(abs(budget)),'maximum_transfer_error',max(abs(transfer_error)), ...
    'maximum_consolidated_budget',max(abs(net_resources)), ...
    'investment_user_price',qI,'consumption_tax_revenue',Tcons, ...
    'intermediate_tax_revenue',Tint,'carbon_tax_revenue',Te, ...
    'investment_tax_revenue',Tx,'government_transfer',transfer, ...
    'government_budget_error',path.transfer-transfer, ...
    'euler_foc',euler,'resource_foc',resource_foc,'household_budget',budget);
metrics=[a.maximum_euler_foc,a.maximum_resource_foc,a.maximum_rent_return_error, ...
    a.maximum_household_budget,a.maximum_transfer_error,a.maximum_consolidated_budget];
a.passed=all(isfinite(metrics))&&max(metrics)<=1e-9;
assert(a.passed,'v2:PolicyEquations','Independent policy FOC/fiscal audit failed');

end
