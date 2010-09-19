function [x res iter]=qve_newton_opt(a,b,tol,maxit)
% [x res iter]=qve_newton_opt(a,M,b,tol,maxit)
% solves a QVE M*x=a+b*(x,x)
% b should be a nxn^2 matrix, b(x,x)=b*kron(x,x)
%
% uses the "traditional" Newton method [Hautphenne, Latouche, Remiche 08]
% "optimized" version --no frills and checks but goes faster
%
% res=residual, iter=number of iterations
% iterates until res<tol or iter>maxit
% (c) f.poloni@sns.it 2009-2010


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
