function X=nr(A,B,C,eps=1.e-12,maxit=inf)
%solves a QBD process X=A+BX+CX^2
%using an inexact Newton method
n=size(A,1);
Ainit=A;
Binit=B;
Cinit=C;
res=inf;
it=0;
X=zeros(n,n);
M=eye(n)-B;
while(res>eps && it<maxit)
  it=it+1;
  W=M\A;
  X=X+W;
  A=C*W*W;
  M=M-2*C*W;
  res=norm(Ainit + Binit*X + Cinit*X*X - X,'fro');
  res,it
end