function M=partialprod(b,x,index);
% function M=partialprod(b,x,index);
% takes b of size n x n^2 and builds b(x,cdot) if index==1 
% or b(cdot,x) if index==2
%
% if index==3, builds M s.t. x'*b(u,v)=v'*M*u

n=size(b,1);
%assert(size(b,2),n^2);
if(index==1)
  M=zeros(n,n);
  for i=1:n
    M=M+b(:,n*(i-1)+(1:n))*x(i);
  end
elseif(index==2)
  M=zeros(n,n);
  for i=1:n
    M=M+b(:,i:n:n^2)*x(i);
  end
elseif(index==3)
  vect=x'*b;
  M=reshape(vect,[n,n]);
else
  error('wrong index!');
end
