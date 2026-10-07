function [loss,audit]=climate_basket_loss(cm,cs,pm,ps,p,lambda)
lc=log([cm*(1-lambda),cs*(1-lambda)]);y=log(pm/ps);
passed=false;iterations=0;
for k=1:60
 [f,d,sm,logE]=residual(y,lc,p);iterations=k;
 assert(isfinite(d)&&d>0&&sm>0&&sm<1,'Climate:PIGLRegularity','Nonregular supporting-price branch');
 if abs(f)<1e-12,passed=true;break;end
 step=-f/d;alpha=1;ok=false;
 for bt=1:30
  yn=y+alpha*step;[fn,dn,sn]=residual(yn,lc,p);
  if isfinite(fn)&&dn>0&&sn>0&&sn<1&&abs(fn)<abs(f),ok=true;break;end
  alpha=alpha/2;
 end
 assert(ok,'Climate:DirectUtility','Supporting-price line search failed');y=yn;
end
assert(passed,'Climate:DirectUtility','Supporting-price solve did not converge');
E=pm*cm+ps*cs;logP=p.omega_m*log(pm)+p.omega_s*log(ps);
Pi=exp(p.zeta_m*log(pm)+p.zeta_s*log(ps));
vb=exp(p.epsilon*(log(E)-logP))/p.epsilon-p.xi*Pi;
vs=exp(p.epsilon*(logE-p.omega_m*y))/p.epsilon-p.xi*exp(p.zeta_m*y);
loss=vb-vs;
assert(isfinite(loss)&&loss>0,'Climate:CE','Scaled basket did not reduce utility');
audit=struct('iterations',iterations,'residual',abs(f),'support_log_price',y,'regularity_derivative',d);
end
function [f,d,sm,logE]=residual(y,lc,p)
v=[y+lc(1),lc(2)];mx=max(v);logE=mx+log(sum(exp(v-mx)));
sm=exp(v(1)-logE);
B=exp((p.zeta_m+p.epsilon*p.omega_m)*y-p.epsilon*logE);
f=sm-p.omega_m-p.xi*p.zeta_m*B;
d=sm*(1-sm)-p.xi*p.zeta_m*B*(p.zeta_m+p.epsilon*p.omega_m-p.epsilon*sm);
end
