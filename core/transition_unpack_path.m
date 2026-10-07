function X = transition_unpack_path(x,q)
assert(isvector(x)&&numel(x)==q.n&&isreal(x)&&all(isfinite(x(:))), ...
    'v2:Coordinates','Invalid coordinate vector');
X=q.reference.*exp(reshape(x,5,q.T+1));
assert(all(isfinite(X(:)))&&all(X(:)>0),'v2:Coordinates','Coordinate overflow/underflow');
end
