function sol=qve_cr(a,M,b);
# solves M*x=a+b*kron(x,x)

n=length(a);

# symmetrizes b
b=(b+b(:,[1 4 7 2 5 8 3 6 9]))/2;

sol=zeros(n,1);
res=1;
a0=a;
M0=M;
while(res>1.E-8)
  iMa=M\a;
  sol=sol+iMa;
  a=b*kron(iMa,iMa);
  for k=1:n
    M=M-2*iMa(k)*b(:,(k-1)*n+(1:n));
  end
  # b does not change
  res=norm(M0*sol-a0-b*kron(sol,sol));
  res
end

end