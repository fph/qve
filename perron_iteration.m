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
% tol=tolerance (e.g. 1.e-15), maxit=max number of iterations (e.g. 10000)
% iterates until res<tol or iter>maxit
% (c) f.poloni@sns.it 2009-2010

n=size(a,1);
e=ones(n,1);
y=a;
x=e-y;res=norm(x-a-b*kron(x,x));

R=partialprod(b,e,1)+partialprod(b,e,2);
[w unused]=perronvector(R',eig_method,tol,e);


norm_helper=w'*(eye(n)-R);
norm_helper2=partialprod(b,w,3);

iter=0;
while(res>n*tol && iter<=maxit)
   iter=iter+1;
   Py=R-partialprod(b,y,1);
   [u lambda]=perronvector(Py,eig_method,tol,y);
   
   t=-(norm_helper*u)/(u'*norm_helper2*u);
   y=t*u;
   
   x=e-y;res=norm(x-a-b*kron(x,x));
   disp(sprintf('%3d residual: %g',iter,res));
end
x=e-y;

