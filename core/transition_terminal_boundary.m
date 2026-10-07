function [f,d] = transition_terminal_boundary(z,s,q)
assert(strcmp(q.boundary_kind,'full_static_abgp_tail'),'v2:Boundary','Unknown closure');
if isfield(q,'investment_policy_kind')&&strcmp(q.investment_policy_kind,'purchase_subsidy')
    assert(all(q.policy.tau_x_path(end-1:end)==0),'v2:Boundary', ...
        'This candidate boundary requires the investment purchase policy to have expired');
end
g=q.growth;tail=(1-g.g_e)*z.R;
f=[log(z.r/g.r_ss);s.capital_residual;(s.e-tail)/tail];
if nargout>1
    d=struct('kind',q.boundary_kind,'residuals',f,'tail_extraction',tail, ...
        'investment_growth_gap',s.I/((g.g_K-1+q.p.delta)*z.K)-1, ...
        'share_corrections',s.shares-[q.p.omega_m,q.p.omega_s], ...
        'feasible',s.feasible,'finite_horizon_validated',false);
end
end
