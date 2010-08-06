function newb=juggleb(b,mode)
%transforms b without modifying the quadratic form t->b(t,t)

n=size(b,1);
if(strcmp(mode,'nothing'))
  newb=b;
elseif(strcmp(mode,'transpose'))
  B=reshape(b,n,n,n);
  B=permute(B,[1 3 2]);
  newb=reshape(B,n,n*n);
elseif(strcmp(mode,'symmetrize'))
  newb=0.5*(b+juggleb(b,'transpose'));
elseif(strcmp(mode,'desymmetrize1'))
  B=reshape(b,n,n,n);
  for j=1:n
    for k=1:j-1
      B(:,j,k)=B(:,j,k)+B(:,k,j);
      B(:,k,j)=0;
    end
  end
  newb=reshape(B,n,n*n);
elseif(strcmp(mode,'desymmetrize2'))
  B=reshape(b,n,n,n);
  for j=1:n
    for k=1:j-1
      B(:,k,j)=B(:,j,k)+B(:,k,j);
      B(:,j,k)=0;
    end
  end
  newb=reshape(B,n,n*n); 
else
  error('Unknown mode');
end
