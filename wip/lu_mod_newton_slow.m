function [x reshistory]=lu_mod_newton_slow(delta,d,q,eps=1.e-12,maxit=inf)
  reshistory=[];
  n=length(delta);
  it=0;
  res=inf;
  P=zeros(n);
  Ptilde=zeros(n);
  for j=1:n
    for i=1:n
      P(i,j)=q(j)/(delta(i)+d(j));
      Ptilde(i,j)=q(j)/(delta(j)+d(i));
    end
  end
  u=ones(n,1);
  v=ones(n,1);
  while(res>eps && it<maxit)
    it=it+1;

    diagRx=1-bdotuv(u,v,P,Ptilde,n);
    diagSx=1./diagRx;
    Sxa=diagSx;
    Jac=eye(2*n)-diag(diagSx)*buvdot(Sxa(1:n),Sxa(n+1:2*n),P,Ptilde,n);
    rhs=[u;v]-Sxa;
    minusdelta=Jac\rhs;
    u=u-minusdelta(1:n);
    v=v-minusdelta(n+1:2*n);
    
    res=norm([u;v].*[P*v;Ptilde*u]+ones(2*n,1)-[u;v],'fro');
    reshistory=[reshistory res];
    res,it
  end
x=[u;v];
end
  
function dd=bdotuv(u,v,P,Ptilde,n)
%returns thex diagonal
  dd=[P*v; Ptilde*u];
end

function M=buvdot(u,v,P,Ptilde,n)
%returns the full (Cauchy-like) matrix
  M=[zeros(n) diag(u)*P; diag(v)*Ptilde zeros(n)];
end