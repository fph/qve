function [x res iter]=perron_newton_eig(a,b,tol,its)
% function x=perron_newton_eig(a,b,tol,its)
%
% as perron_newton, but 
%
% (c) f.poloni@sns.it 2009-2010

n=size(a,1);
e=ones(n,1);
iter=0;

R=partialprod(b,e,1)+partialprod(b,e,2);
[V Lambda]=eig(R');
[useless,j]=max(real(diag(Lambda)));
w=V(:,j);

norm_helper=w'*(eye(n)-R);
norm_helper2=partialprod(b,w,3);

v=e; %"seed" for perronvector for the left eigenvector
y=a;

for iter=1:its
   P=R-partialprod(b,y,1);
   [V Lambda]=eig(P);
   [useless,j]=max(real(diag(Lambda)));
   u=V(:,j);%u=real(u/u(1));
   lambda=Lambda(j,j);
   v=zeros(1,n);v(1,j)=1;v=v/V;v=v';v=real(v/v(1)); %instead of doing a second eig(), extracts the needed left eigenvector v by inverting the the right-eigenvector matrix
   
   t=-(norm_helper*u)/(u'*norm_helper2*u);
   
   u=t*u;
   Temp=pinv(P-lambda*eye(n))*(eye(n)-u*v'/(v'*u))*partialprod(b,u,2);
   sigmaT=norm_helper+u'*(norm_helper2+norm_helper2');
   Jac=(eye(n)-u*sigmaT/(sigmaT*u))*Temp;
   y=y-(eye(n)-Jac)\(y-u);

   x=e-y;res=norm(x-a-b*kron(x,x));
   disp(sprintf('%3d residual: %g',iter,res));
   if(res<n*tol)
       break;
   end
end
x=e-y;

