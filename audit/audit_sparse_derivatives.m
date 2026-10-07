function a = audit_sparse_derivatives(x,q)
x=x+0.001*sin((1:q.n)');
[J,info]=transition_sparse_jacobian(x,q,2e-5);
Jhalf=transition_sparse_jacobian(x,q,1e-5);
a=struct('passed',false,'step_error',norm(J-Jhalf,inf)/max(1,norm(Jhalf,inf)), ...
    'column_error',0,'off_pattern_error',0,'direction_error',0,'info',info);
last=q.policy_specification.start+q.policy_specification.duration-1;
nodes=unique([0,max(0,last-1),last,last+1,floor(q.T/2),q.T]);nodes=nodes(nodes<=q.T);cols=[];
for t=nodes,cols=[cols,5*t+(1:5)];end 
for c=cols
    h=8e-6*max(1,abs(x(c)));xp=x;xm=x;xp(c)=xp(c)+h;xm(c)=xm(c)-h;
    v=(transition_global_residual(xp,q)-transition_global_residual(xm,q))/(2*h);
    a.column_error=max(a.column_error,norm(v-full(J(:,c)),inf)/max(1,norm(v,inf)));
    a.off_pattern_error=max(a.off_pattern_error,norm(v(~logical(q.pattern(:,c))),inf));
end
for k=1:3
    v=sin((1:q.n)'*(0.37+k*0.29));v=v/norm(v,inf);h=8e-6;
    obs=(transition_global_residual(x+h*v,q)-transition_global_residual(x-h*v,q))/(2*h);
    a.direction_error=max(a.direction_error,norm(obs-J*v,inf)/max(1,norm(obs,inf)));
end
a.selected_columns=cols;
a.passed=issparse(J)&&a.step_error<2e-7&&a.column_error<2e-7 ...
    &&a.direction_error<2e-7&&a.off_pattern_error<1e-10&&sprank(q.pattern)==q.n;
assert(a.passed,'v2:Derivative','Sparse policy derivative audit failed');

end
