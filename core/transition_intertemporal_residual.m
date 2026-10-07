function [f,d] = transition_intertemporal_residual(E,En,P,Pn,rn,h,hn,p,tx,txn,mode,purchase)
mu=(p.epsilon-1)*log(E)-p.epsilon*log(P);
mun=(p.epsilon-1)*log(En)-p.epsilon*log(Pn);
gross=1-p.delta+rn;
if purchase
    assert(strcmp(mode,'paper'),'v2:Scope','Purchase subsidy requires paper equations');
    qI=1+tx;qIn=1+txn;payoff=rn+(1-p.delta)*qIn;
    assert(qI>0&&qIn>0&&payoff>0,'v2:Return','Invalid purchase-policy return');
    gross=payoff/qI;
elseif strcmp(mode,'legacy')
    mu=mu-log1p(tx);mun=mun-log1p(txn);
else
    assert(tx==0&&txn==0,'v2:Scope','Unspecified paper investment wedge');
end
assert(gross>0&&isfinite(gross),'v2:Return','Invalid gross return');
f=[mu-mun-log(p.beta)-log(gross);log(hn)-log(h)-log(gross)];
if nargout>1
    d=struct('log_mu',mu,'log_mu_next',mun,'capital_gross_return',gross, ...
        'purchase_policy',purchase);
end
end
