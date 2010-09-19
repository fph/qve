function r=convradius(a,b)
% function convradius(a,b)
% return the value of rho(Jacobian(minimal sol'n))

n=length(a);
e=ones(n,1);
x=perron_iteration_newnorm_eig(a,b,1.e-14,9999999);
M=partialprod(b,e,1)+partialprod(b,e,2);
w=perronvector(M','eig',1.e-14);
Jac=perron_jacobian(a,b,'eig',ones(size(a))-x,w);
r=max(abs(eig(Jac)));
