function [x res iter]=perron_newton(a,b,eig_method,tol,its)
%
% function x=perron_newton(a,b,eig_method,tol,its)
%
% computes the minimal solution of a MBT using the Perron Newton method
% tol=max. relative residual allowed
% its=max. iterations (may be inf)

n=size(a,1);
e=ones(n,1);
y=a;
x=e-y;res=norm(x-a-b*kron(x,x));
iter=0;
%disp(sprintf('%3d residual: %g',iter,res));

v=rand(n,1); %starting vector for the left eigvec
for iter=1:its
   P=partialprod(b,e-y,1)+partialprod(b,e,2);
   [u lambda]=perronvector(P,eig_method,tol,y);
   [v lambda2]=perronvector(P',eig_method,tol,v);

   v1=u-b*kron(e,u)-b*kron(u,e);
   v2=b*kron(u,u);
   t=-(e'*v1)/(e'*v2);
   u=t*u;
   R=pinv(P-lambda*eye(n))*(eye(n)-u*v'/(v'*u))*partialprod(b,u,2);
   sigmaT=e'*(eye(n)-partialprod(b,e-u,1)-partialprod(b,e-u,2));
   Jac=(eye(n)-u*sigmaT/(sigmaT*u))*R;

   y=y-(eye(n)-Jac)\(y-u);

   x=e-y;res=norm(x-a-b*kron(x,x));
   %disp(sprintf('%3d residual: %g',iter,res));
   if(res<n*tol)
       break;
   end
end
x=e-y;
