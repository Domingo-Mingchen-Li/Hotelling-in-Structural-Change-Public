function path = transition_export_path(x,q)
X=transition_unpack_path(x,q);path=struct();
for k=1:5,path.(q.variable_names{k})=X(k,:);end
for k=1:4
    names={'A_m','A_s','A_x','A_e'};path.(names{k})=q.technology(k,:);
end
for j=1:q.T+1
    z=struct('K',X(1,j),'R',X(2,j),'r',X(3,j),'h',X(4,j),'E',X(5,j));
    tech=struct('A_m',q.technology(1,j),'A_s',q.technology(2,j), ...
        'A_x',q.technology(3,j),'A_e',q.technology(4,j));
    pol=struct('tau_c_m',q.policy.tau_c_m_path(j),'tau_c_s',q.policy.tau_c_s_path(j), ...
        'tau_int',q.policy.tau_int_path(j),'tau_e',q.policy.tau_e_path(j),'tau_x',q.policy.tau_x_path(j));
    if isfield(q,'investment_policy_kind'),pol.investment_policy_kind=q.investment_policy_kind;end
    s=transition_static_kernel(z,tech,pol,q.p,q.mode);
    for name={'w','p_m','p_s','p_int','p_int_user','P_star','I','c_m','c_s','K_d','transfer','investment_user_price','investment_transfer','budget_raw'}
        path.(name{1})(j)=s.(name{1});
    end
    path.R_input(j)=s.e;path.s_m(j)=s.shares(1);path.s_s(j)=s.shares(2);
    for k=1:4
        sectors={'m','s','x','e'};
        path.(['k_' sectors{k}])(j)=s.k(k);path.(['l_' sectors{k}])(j)=s.l(k);
    end
end
end
