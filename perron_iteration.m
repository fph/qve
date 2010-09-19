function [x res iter]=perron_iteration(a,b,eig_method,tol,maxit)
%
% [x res iter]=perron_iteration(a,b,eig_method,tol,maxit)
% solves a QVE M*x=a+b*(x,x)
% b should be a nxn^2 matrix, b(x,x)=b*kron(x,x)
% assumes M=I and that x=e is a (nonnecessarily minimal) solution
%
% uses the Perron iteration [Meini, Poloni '10 arXiv]
%
% eig_method is as in perronvector.m
%
% res=residual, iter=number of iterations
% iterates until res<tol or iter>maxit
% (c) f.poloni@sns.it 2009-2010

if(not(exist('maxit','var')))
   maxit=inf;
end

n=size(a,1);
e=ones(n,1);
y=a;
x=e-y;res=norm(x-a-b*kron(x,x));
iter=0;
%disp(sprintf('%3d residual: %g',iter,res));

R=partialprod(b,e,1)+partialprod(b,e,2);

iter=0;
while(res>n*tol && iter<=maxit)
   iter=iter+1;
   Py=R-partialprod(b,y,1);
   [u lambda]=perronvector(Py,eig_method,tol,y);
   v1=u-b*kron(e,u)-b*kron(u,e);
   v2=b*kron(u,u);
   t=-(e'*v1)/(e'*v2);
   y=t*u;
   
   x=e-y;res=norm(x-a-b*kron(x,x));
   disp(sprintf('%3d residual: %g',iter,res));
end
x=e-y;
