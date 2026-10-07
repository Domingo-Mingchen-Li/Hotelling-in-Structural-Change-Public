function warm = select_initial_path(pool,model,id,T)
idx=find(cellfun(@(s)strcmp(s.model_mode,model)&&strcmp(s.case_name,id),pool));
if isempty(idx),idx=find(cellfun(@(s)strcmp(s.model_mode,model)&&strcmp(s.case_name,'baseline'),pool));end
assert(~isempty(idx),'v2:WarmStart','No %s initial path available',model);
h=cellfun(@(s)s.T,pool(idx));long=find(h>=T);
if isempty(long),[~,k]=max(h);else,[~,kk]=min(h(long));k=long(kk);end
warm=pool{idx(k)};
assert(warm.validation.passed&&strcmp(warm.boundary_kind,'full_static_abgp_tail'), ...
    'v2:WarmStart','Initial path lacks finite-horizon validation');
end
