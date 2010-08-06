function Y=ginv(A)
%group inverse, as in Meyer-Stewart
n=size(A,1);
U=orth(A);
X=null(A);
r=size(U,2);
P=[X U];
C=U'*A*U;
bigC=[zeros(n-r,n-r) zeros(n-r,r);zeros(r,n-r) C];
Y=P*bigC/P;