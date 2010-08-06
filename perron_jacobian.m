function [Jac dyad p1 dyad2 p2 sigmaT pp1 v]=perron_jacobian(a,b,eig_method,y,w)
%computes the Jacobian of the Perron iteration in a point y
%XX: test only
   n=size(a,1);
   e=ones(n,1);
   tol=1.e-13;
   
   v=rand(n,1);
   P=partialprod(b,e-y,1)+partialprod(b,e,2);
   [u lambda]=perronvector(P,eig_method,tol,y);
   [v lambda2]=perronvector(P',eig_method,tol,v);
   
   v1=u-b*kron(e,u)-b*kron(u,e);
   v2=b*kron(u,u);
   t=-(w'*v1)/(w'*v2);
   u=t*u;
   R=pinv(P-lambda*eye(n))*(eye(n)-u*v'/(v'*u))*partialprod(b,u,2);
   sigmaT=w'*(eye(n)-partialprod(b,e-u,1)-partialprod(b,e-u,2));
   Jac=(eye(n)-u*sigmaT/(sigmaT*u))*R;

   dyad=eye(n)-u*sigmaT/(sigmaT*u);
   p1=pinv(P-lambda*eye(n));
   pp1=P-lambda*eye(n);
   dyad2=(eye(n)-u*v'/(v'*u));
   p2=partialprod(b,u,2);
end