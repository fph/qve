function [x res iter]=perron_newton_newnorm_eig(a,b,tol,its)
%
% function x=perron_newton(a,b,tol,its)
%
% computes the minimal solution of a MBT using the Perron Newton method
% tol=max. relative residual allowed
% its=max. iterations (may be inf)

%uses eig to compute ev.
%normalizes orthogonalizing wrt the left eigv of M

n=size(a,1);
e=ones(n,1);
y=a;
x=e-y;res=norm(x-a-b*kron(x,x));
iter=0;
%disp(sprintf('%3d residual: %g',iter,res));

M=partialprod(b,e,1)+partialprod(b,e,2);
[V Lambda]=eig(M');
[useless,j]=max(real(diag(Lambda)));
w=V(:,j);

v=rand(n,1); %starting vector for the left eigvec
for iter=1:its
   P=partialprod(b,e-y,1)+partialprod(b,e,2);
   
   [V Lambda]=eig(P);
   [useless,j]=max(real(diag(Lambda)));
   u=V(:,j);u=real(u/u(1));
   
   lambda=Lambda(j,j);
   v=zeros(1,n);v(1,j)=1;v=v/V;v=v';v=real(v/v(1));
   
   %[vtrue lambda2]=perronvector(P','eig',tol,v);
   %v./vtrue

   v1=u-b*kron(e,u)-b*kron(u,e);
   v2=b*kron(u,u);
   t=-(w'*v1)/(w'*v2);
   u=t*u;
   R=pinv(P-lambda*eye(n))*(eye(n)-u*v'/(v'*u))*partialprod(b,u,2);
   sigmaT=w'*(eye(n)-partialprod(b,e-u,1)-partialprod(b,e-u,2));
   Jac=(eye(n)-u*sigmaT/(sigmaT*u))*R;

   y=y-(eye(n)-Jac)\(y-u);

   x=e-y;res=norm(x-a-b*kron(x,x));
   disp(sprintf('%3d residual: %g',iter,res));
   if(res<n*tol)
       break;
   end
end
x=e-y;
