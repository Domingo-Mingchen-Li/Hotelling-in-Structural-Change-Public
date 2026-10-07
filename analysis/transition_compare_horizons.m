function c = transition_compare_horizons(a,b,window_end)
assert(window_end<=min(a.T,b.T),'v2:Comparison','Window outside common horizon');
idx=1:window_end+1;
names={'K','R','r','h','E','I','R_input','c_m','c_s','w','p_m','p_s','P_star'};
errors=zeros(size(names));
for k=1:numel(names)
    key=names{k};ya=a.path.(key)(idx);yb=b.path.(key)(idx);
    assert(all(yb>0),'v2:Comparison','Nonpositive comparison level');
    errors(k)=max(abs(ya./yb-1));
end
share_error=max(abs(a.path.s_m(idx)-b.path.s_m(idx)));
c=struct('shorter_T',a.T,'longer_T',b.T,'window_end',window_end, ...
    'variable_names',{names},'maximum_relative_differences',errors, ...
    'maximum_over_variables',max(errors),'maximum_share_difference',share_error, ...
    'share_difference_percentage_points',100*share_error);
end
