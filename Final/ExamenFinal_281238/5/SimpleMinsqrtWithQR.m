A=[3,-1;4,2;0,1];
b=[0;2;1];
[Q,R]=qr(A);
b1=Q'*b;;
sol=R\b1