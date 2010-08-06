function x=order(a,M,b,epsilon=1.e-8)
n=length(a);
x=zeros(n,1);
xold=37;
it=0;
while(norm(x-xold)>epsilon);
  it=it+1;
  xold=x;
  N=M-partialprod(b,xold,1);
  x=N\a;
  it,norm(x-xold)
end