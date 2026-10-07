function [P,colors] = transition_jacobian_structure(T)
n=5*(T+1);P=sparse(n,n);
for t=0:T-1
    c=5*t;next=c+5;row=5*t;
    P(row+1,c+[1,3,4,5])=1;                 
    P(row+2,[c+[1,3,4,5],next+1])=1;        
    P(row+3,[c+[2,3,4,5],next+2])=1;        
    P(row+4,[c+[3,4,5],next+[3,4,5]])=1;    
    P(row+5,[c+4,next+[3,4]])=1;            
end
row=5*T;P(row+1,1)=1;P(row+2,2)=1;
P(row+3,5*T+3)=1;
P(row+4,5*T+[1,3,4,5])=1;
P(row+5,5*T+[2,3,4,5])=1;
P=spones(P);
colors=zeros(n,1);
for c=1:n
    rows=find(P(:,c));neighbors=find(any(P(rows,1:c-1),1));
    used=colors(neighbors);color=1;
    while any(used==color),color=color+1;end
    colors(c)=color;
end
for color=1:max(colors)
    assert(all(full(sum(P(:,colors==color),2))<=1),'v2:Color','Overlapping color columns');
end
end
