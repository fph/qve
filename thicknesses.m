function [x res iter]=thicknesses(a,M,b,tol,maxit)

if(not(exist('maxit','var')))
   maxit=inf;
end

n=length(a);

x=zeros(n,1);
res=inf;
iter=0;
while(res>n*tol && iter<=maxit)
  iter=iter+1;
  mat=M-partialprod(b,x,mod(iter,2)+1);
  x=mat\a;
  res=norm(M*x-a-b*kron(x,x));
%  disp(sprintf('%3d residual: %g',iter,res));
end



