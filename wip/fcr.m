function X=fcr(A,B,C,eps=1.e-12,maxit=inf)
%solves a QBD process X=A+BX+CX^2
%new "greedy" LR

n=size(A,1);

Ainit=A;
Binit=B;
Cinit=C;

A=(eye(n)-B)\A;
C=(eye(n)-B)\C;

X=A;

PC=C;

res=inf;
it=0;
while(res>eps && it<maxit)
  it=it+1;
  %R should be <=X (closer possible) and commute with X
  %R<=X iff R<=A+CR^2
  %we shall look for R=alpha*I
  v=1-4*diag(A).*diag(C);
  A
  C
  alpha=zeros(n,1);
  for i=1:n
    if(v(i)>=0)
      alpha(i)=( 1-sqrt(v(i)) ) /2;
    else
      alpha(i)=inf;
    end
  end
  alpha=max(alpha);
  alpha2=alpha^2;
  %first correction to X
  X=X+PC*alpha2;
  %creates the residual equation for X2-R2
  Atilde=A+C*alpha2;
  M=eye(n)-Atilde*C-C*Atilde;
  A=M\ (Atilde^2-alpha2*eye(n));
  C=M\ (C^2);
  X=X+PC*A;
  PC=PC*C;

  res=norm(Ainit + Binit*X + Cinit*X*X - X,'fro');
  res,it
end
