function [value,detail]=climate_environment_tail(C,H,qT,ge,beta,m,Mref,eta)
A=[diag(m.rho),m.a;zeros(1,4),ge];b=[ones(4,1);0];
h=[H(:);0];f=[C(:)-H(:);qT];
assert(Mref>0&&ge>0&&ge<1&&beta>0&&beta<1,'Climate:Tail','Invalid tail');
if eta==1
 P=[];v=(eye(5)-beta*A')\b;
 value=v'*f/Mref;
 gradient=v/Mref;
else
 P=reshape((eye(25)-beta*kron(A',A'))\reshape(b*b',25,1),5,5);
 P=(P+P')/2;
 value=(f'*P*f+2*h'*P*f)/Mref^2;
 gradient=2*P*(h+f)/Mref^2;
 v=[];
end
assert(isfinite(value)&&value>=-1e-10,'Climate:Tail','Invalid damage continuation');
detail=struct('A',A,'P',P,'linear_value',v,'gradient_total_state',gradient, ...
 'terminal_future_state',f,'terminal_inherited_state',h,'exact_given_geometric_tail',true);
end
