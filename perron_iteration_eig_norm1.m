function [x res iter]=perron_iteration_eig(a,b,method,tol,maxit)
%
% function x=perron_iteration_eig(a,b,method,tol,its)
%
% computes the minimal solution of a MBT using the Perron iteration
% tol=max. relative residual allowed
% its=max. iterations (may be inf)
% uses "eig" for eigenvalues
%normalizes with left Perron vector of R

if(not(exist('maxit','var')))
   maxit=inf;
end

n=size(a,1);
e=ones(n,1);
y=a;
x=e-y;res=norm(x-a-b*kron(x,x));

R=partialprod(b,e,1)+partialprod(b,e,2);
[V Lambda]=eig(R');
[useless,j]=max(real(diag(Lambda)));
w=V(:,j);


norm_helper=w'*(eye(n)-R);

iter=0;
while(res>n*tol && iter<=maxit)
   iter=iter+1;
   Py=R-partialprod(b,y,1);
   [V Lambda]=eig(Py);
   [useless,j]=max(real(diag(Lambda)));
   u=V(:,j);
   
   v2=b*kron(u,u);
   t=-(norm_helper*u)/(w'*v2);
   y=t*u;
   
   x=e-y;res=norm(x-a-b*kron(x,x));
   disp(sprintf('%3d residual: %g',iter,res));
end
x=e-y;