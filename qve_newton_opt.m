function [x res iter]=qve_newton_opt(a,b,tol,maxit)
n=length(a);

x=zeros(n,1);
res=inf;
iter=0;
bxx=b*kron(x,x);
while(res>n*tol && iter<=maxit)
  iter=iter+1;
  rhs=a-bxx;
  mat=eye(n)-partialprod(b,x,1)-partialprod(b,x,2);
  x=mat\rhs;
  bxx=b*kron(x,x);
  res=norm(x-a-bxx);
%  disp(sprintf('%3d residual: %g',iter,res));
end
