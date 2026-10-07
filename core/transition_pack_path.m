function x = transition_pack_path(X,q)
assert(isequal(size(X),[5,q.T+1])&&isreal(X)&&all(isfinite(X(:))) ...
    &&all(X(:)>0),'v2:Path','Expected positive finite 5-by-(T+1) path');
U=log(X./q.reference);x=U(:);
end
