function [x res iter]=qve_newton(a,M,b,tol,maxit)

if(not(exist('maxit','var')))
   maxit=inf;
end

n=length(a);

x=zeros(n,1);
res=inf;
iter=0;
oldx=inf;
while(res>n*tol && iter<=maxit)
  oldx=x;
  iter=iter+1;
  rhs=a-b*kron(x,x);
  mat=M-partialprod(b,x,1)-partialprod(b,x,2);
  x=mat\rhs;
  res=norm(M*x-a-b*kron(x,x));
  disp(sprintf('%3d residual: %g',iter,res));
end
