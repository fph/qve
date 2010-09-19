function X=lr(A,B,C,eps=1.e-12,maxit=inf)
%solves a QBD process X=A+BX+CX^2
%logarithmic reduction

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

  %creates the residual equation
  M=eye(n)-A*C-C*A;
  A=M\ (A^2);
  C=M\ (C^2);
  X=X+PC*A;
  PC=PC*C;

  res=norm(Ainit + Binit*X + Cinit*X*X - X,'fro');
  res,it
end