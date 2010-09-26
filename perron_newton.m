function [x res iter]=perron_newton(a,b,eig_method,tol,its)
%
% function x=perron_newton(a,b,eig_method,tol,its)
%
% solves a QVE M*x=a+b*(x,x)
% b should be a nxn^2 matrix, b(x,x)=b*kron(x,x)
% assumes M=I and that x=e is a (nonnecessarily minimal) solution
%
% uses the Perron-Newton method [Meini, Poloni '10 arXiv]
%
% eig_method is as in perronvector.m
%
% res=residual, iter=number of iterations
% tol=tolerance (e.g. 1.e-15), maxit=max number of iterations (e.g. 10000)
% iterates until res<tol or iter>maxit
% (c) f.poloni@sns.it 2009-2010

n=size(a,1);
e=ones(n,1);
iter=0;

R=partialprod(b,e,1)+partialprod(b,e,2);
[w unused]=perronvector(R',eig_method,tol,e);

norm_helper=w'*(eye(n)-R);
norm_helper2=partialprod(b,w,3);

v=e; %"seed" for perronvector for the left eigenvector
y=a;

for iter=1:its
   P=R-partialprod(b,y,1);
   [u lambda]=perronvector(P,eig_method,tol,y);
   [v lambda2]=perronvector(P',eig_method,tol,v);
   
   t=-(norm_helper*u)/(u'*norm_helper2*u);
   
   u=t*u;
   
   Temp=pinv(P-lambda*eye(n));
   Temp=Temp-u*(v'*Temp)/(v'*u);
   Temp=Temp*partialprod(b,u,2);
   %Temp=pinv(P-lambda*eye(n))*(eye(n)-u*v'/(v'*u))*partialprod(b,u,2);
   sigmaT=norm_helper+u'*(norm_helper2+norm_helper2');
   %Jac=(eye(n)-u*sigmaT/(sigmaT*u))*Temp;
   Jac=Temp-u*(sigmaT*Temp)/(sigmaT*u);
   y=y-(eye(n)-Jac)\(y-u);

   x=e-y;res=norm(x-a-b*kron(x,x));
   disp(sprintf('%3d residual: %g',iter,res));
   if(res<n*tol)
       break;
   end
end
x=e-y;

