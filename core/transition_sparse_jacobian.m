function [J,info] = transition_sparse_jacobian(x,q,relative_step)
if nargin<3,relative_step=2e-5;end
assert(isscalar(relative_step)&&isfinite(relative_step)&&relative_step>0, ...
    'v2:Derivative','Invalid finite-difference step');
x=x(:);steps=relative_step*max(1,abs(x));
[rows,cols]=find(q.pattern);values=zeros(size(rows));
for color=1:max(q.colors)
    members=find(q.colors==color);xp=x;xm=x;
    xp(members)=xp(members)+steps(members);xm(members)=xm(members)-steps(members);
    fp=transition_global_residual(xp,q);fm=transition_global_residual(xm,q);
    selected=find(q.colors(cols)==color);
    values(selected)=(fp(rows(selected))-fm(rows(selected)))./(2*steps(cols(selected)));
end
assert(all(isfinite(values)),'v2:Derivative','Nonfinite Jacobian entries');
J=sparse(rows,cols,values,q.n,q.n);
info=struct('method','colored central differences','relative_step',relative_step, ...
    'colors',max(q.colors),'residual_evaluations',2*max(q.colors), ...
    'structural_nnz',nnz(q.pattern),'numeric_nnz',nnz(J));
end
