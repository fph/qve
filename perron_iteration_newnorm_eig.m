function [x res iter]=perron_iteration_newnorm_eig(a,b,tol,its)
%
% function x=perron_iteration(a,b,tol,its)
%
% computes the minimal solution of a MBT using the Perron iteration
% tol=max. relative residual allowed
% its=max. iterations (may be inf)
% uses "eig" for eigenvalues
% normalizes using the left eigv of M

n=size(a,1);
e=ones(n,1);

M=partialprod(b,e,1)+partialprod(b,e,2);
[V Lambda]=eig(M');
[useless,j]=max(real(diag(Lambda)));
w=V(:,j);

y=a;
x=e-y;res=norm(x-a-b*kron(x,x));
iter=0;
%disp(sprintf('%3d residual: %g',iter,res));

for iter=1:its
   Py=partialprod(b,e-y,1)+partialprod(b,e,2);
   [V Lambda]=eig(Py);
   [useless,j]=max(real(diag(Lambda)));
   u=V(:,j);
   
   v1=u-b*kron(e,u)-b*kron(u,e);
   v2=b*kron(u,u);
   t=-(w'*v1)/(w'*v2);
   y=t*u;
   
   x=e-y;res=norm(x-a-b*kron(x,x));
%  disp(sprintf('%3d residual: %g',iter,res));
   if(res<n*tol)
       break;
   end
end
x=e-y;