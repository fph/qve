function [x res iter]=perron_iteration(a,b,eig_method,tol,its)
%
% function x=perron_iteration(a,b,eig_method,tol,its)
%
% computes the minimal solution of a MBT using the Perron iteration
% tol=max. relative residual allowed
% its=max. iterations (may be inf)

n=size(a,1);
e=ones(n,1);
y=a;
x=e-y;res=norm(x-a-b*kron(x,x));
iter=0;
%disp(sprintf('%3d residual: %g',iter,res));

for iter=1:its
   Py=partialprod(b,e-y,1)+partialprod(b,e,2);
   [u lambda]=perronvector(Py,eig_method,tol,y);
   v1=u-b*kron(e,u)-b*kron(u,e);
   v2=b*kron(u,u);
   t=-(e'*v1)/(e'*v2);
   y=t*u;
   
   x=e-y;res=norm(x-a-b*kron(x,x));
%   disp(sprintf('%3d residual: %g',iter,res));
   if(res<n*tol)
       break;
   end
end
x=e-y;
